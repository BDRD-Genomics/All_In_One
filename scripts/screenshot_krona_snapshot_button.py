#!/usr/bin/env python3

import argparse
import re
from pathlib import Path
from urllib.parse import unquote

import cairosvg
from playwright.sync_api import sync_playwright, TimeoutError as PlaywrightTimeoutError


def add_white_background(svg_text: str) -> str:
    """
    Add a white background rectangle as the first child of the SVG.
    This prevents transparent backgrounds when converting to PNG or inserting into DOCX.
    """
    if 'id="white_background"' in svg_text:
        return svg_text

    return re.sub(
        r'(<svg\b[^>]*>)',
        r'\1\n<rect id="white_background" width="100%" height="100%" fill="white"/>',
        svg_text,
        count=1,
    )


def extract_svg_from_data_uri(data_uri: str) -> str | None:
    """
    Extract SVG XML from a data:image/svg+xml URI.
    """
    if not data_uri or not data_uri.startswith("data:image/svg+xml"):
        return None

    comma_index = data_uri.find(",")
    if comma_index < 0:
        return None

    return unquote(data_uri[comma_index + 1:])


def scale_svg_font_sizes(
    svg_text: str,
    scale: float | None = None,
    fixed_px: float | None = None,
) -> str:
    """
    Scale or set font sizes in an SVG.

    Examples:
      --font-scale 1.5  increases existing font sizes by 50%
      --font-size 18    sets all detected font sizes to 18px

    If neither scale nor fixed_px is provided, the SVG is returned unchanged.
    """
    if scale is None and fixed_px is None:
        return svg_text

    if scale is not None and scale <= 0:
        raise ValueError("--font-scale must be greater than 0")

    if fixed_px is not None and fixed_px <= 0:
        raise ValueError("--font-size must be greater than 0")

    def replace_style_font_size(match):
        value = float(match.group(1))
        unit = match.group(2) or "px"

        if fixed_px is not None:
            return f"font-size:{fixed_px:g}px"

        return f"font-size:{value * scale:g}{unit}"

    # Handles inline CSS styles like:
    # font-size:12px
    # font-size: 12
    # font-size:10pt
    svg_text = re.sub(
        r"font-size\s*:\s*([0-9.]+)\s*(px|pt|em|rem)?",
        replace_style_font_size,
        svg_text,
        flags=re.IGNORECASE,
    )

    def replace_attr_font_size(match):
        value = float(match.group(1))
        unit = match.group(2) or "px"

        if fixed_px is not None:
            return f'font-size="{fixed_px:g}px"'

        return f'font-size="{value * scale:g}{unit}"'

    # Handles XML attributes like:
    # font-size="12px"
    # font-size="12"
    # font-size="10pt"
    svg_text = re.sub(
        r'font-size="([0-9.]+)\s*(px|pt|em|rem)?"',
        replace_attr_font_size,
        svg_text,
        flags=re.IGNORECASE,
    )

    return svg_text


def svg_to_png_with_cairosvg(
    svg_path: Path,
    png_path: Path,
    output_width: int | None = None,
    output_height: int | None = None,
):
    """
    Convert SVG to PNG using CairoSVG.
    """
    png_path.parent.mkdir(parents=True, exist_ok=True)

    kwargs = {
        "url": str(svg_path),
        "write_to": str(png_path),
        "background_color": "white",
    }

    if output_width:
        kwargs["output_width"] = output_width

    if output_height:
        kwargs["output_height"] = output_height

    cairosvg.svg2png(**kwargs)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Open a Krona HTML file, click Snapshot, wait for the popup, "
            "click Download Snapshot, save SVG, optionally adjust font size, "
            "and optionally convert SVG to PNG using CairoSVG."
        )
    )

    parser.add_argument("-i", "--input", required=True, help="Input Krona HTML file")
    parser.add_argument("-o", "--output", required=True, help="Output SVG file")

    parser.add_argument(
        "-p",
        "--png-output",
        default=None,
        help="Optional output PNG file converted from the saved SVG using CairoSVG",
    )

    # Keep -j as a compatibility alias, but treat it as PNG output.
    # This lets old commands fail less abruptly if you update them gradually.
    parser.add_argument(
        "-j",
        "--jpeg-output",
        default=None,
        help=(
            "Deprecated: use --png-output. "
            "If provided, this path will be used as PNG output if it ends in .png."
        ),
    )

    parser.add_argument(
        "--width",
        type=int,
        default=1600,
        help="Browser viewport width for Krona rendering",
    )

    parser.add_argument(
        "--height",
        type=int,
        default=1200,
        help="Browser viewport height for Krona rendering",
    )

    parser.add_argument(
        "--wait-ms",
        type=int,
        default=8000,
        help="Wait after Krona page load",
    )

    parser.add_argument(
        "--popup-timeout-ms",
        type=int,
        default=30000,
        help="Timeout waiting for the Krona snapshot popup",
    )

    parser.add_argument(
        "--download-timeout-ms",
        type=int,
        default=30000,
        help="Timeout waiting for the SVG download",
    )

    parser.add_argument(
        "--png-width",
        type=int,
        default=2200,
        help="Output PNG width for CairoSVG conversion",
    )

    parser.add_argument(
        "--png-height",
        type=int,
        default=None,
        help="Optional output PNG height for CairoSVG conversion",
    )

    parser.add_argument(
        "--font-scale",
        type=float,
        default=None,
        help=(
            "Scale detected SVG font sizes. "
            "Example: --font-scale 1.5 makes fonts 50 percent larger."
        ),
    )

    parser.add_argument(
        "--font-size",
        type=float,
        default=None,
        help=(
            "Set detected SVG font sizes to a fixed pixel size. "
            "Example: --font-size 18"
        ),
    )

    parser.add_argument(
        "--no-white-background",
        action="store_true",
        help="Do not inject a white background rectangle into the SVG",
    )

    args = parser.parse_args()

    html_path = Path(args.input).resolve()
    out_path = Path(args.output).resolve()

    png_output = args.png_output
    if not png_output and args.jpeg_output:
        png_output = args.jpeg_output

    png_path = Path(png_output).resolve() if png_output else None

    if not html_path.exists():
        raise FileNotFoundError(f"Input HTML does not exist: {html_path}")

    if out_path.suffix.lower() != ".svg":
        raise ValueError(f"Output file should end in .svg: {out_path}")

    if png_path and png_path.suffix.lower() != ".png":
        raise ValueError(
            f"PNG output should end in .png: {png_path}. "
            "Use -p sample.krona.snapshot.png instead of -j sample.jpeg."
        )

    if args.font_scale is not None and args.font_size is not None:
        raise ValueError(
            "Use either --font-scale or --font-size, not both. "
            "Recommended: start with --font-scale 1.4"
        )

    out_path.parent.mkdir(parents=True, exist_ok=True)

    svg_text = None
    triggered = "not_triggered"
    clicked = "not_clicked"
    suggested_filename = None

    with sync_playwright() as p:
        browser = p.chromium.launch(
            headless=True,
            args=[
                "--no-sandbox",
                "--disable-dev-shm-usage",
                "--disable-gpu",
                "--allow-file-access-from-files",
            ],
        )

        context = browser.new_context(
            viewport={"width": args.width, "height": args.height},
            device_scale_factor=1,
            accept_downloads=True,
        )

        page = context.new_page()

        page.goto(html_path.as_uri(), wait_until="load")
        page.wait_for_timeout(args.wait_ms)

        try:
            with page.expect_popup(timeout=args.popup_timeout_ms) as popup_info:
                triggered = page.evaluate(
                    """
                    () => {
                        if (typeof snapshot === 'function') {
                            snapshot();
                            return 'snapshot_function';
                        }

                        const snapshotInput = document.querySelector('#snapshot');
                        if (snapshotInput) {
                            snapshotInput.click();
                            return 'snapshot_button_id';
                        }

                        if (typeof snapshotButton !== 'undefined' && snapshotButton) {
                            snapshotButton.click();
                            return 'snapshotButton_variable';
                        }

                        const candidates = Array.from(document.querySelectorAll('input, button, a, img'));
                        const snap = candidates.find(el => {
                            const text = (
                                el.innerText ||
                                el.value ||
                                el.title ||
                                el.alt ||
                                el.getAttribute('aria-label') ||
                                ''
                            ).toLowerCase();

                            return text.includes('snapshot');
                        });

                        if (snap) {
                            snap.click();
                            return 'dom_snapshot_click';
                        }

                        return 'not_found';
                    }
                    """
                )

            popup = popup_info.value
            popup.wait_for_load_state("load", timeout=args.popup_timeout_ms)
            popup.wait_for_timeout(1000)

        except PlaywrightTimeoutError:
            browser.close()
            raise RuntimeError(
                "Snapshot was triggered, but no popup window appeared."
            )

        try:
            with popup.expect_download(timeout=args.download_timeout_ms) as download_info:
                clicked = popup.evaluate(
                    """
                    () => {
                        const links = Array.from(document.querySelectorAll('a'));
                        const downloadLink = links.find(a => {
                            const text = (a.innerText || a.textContent || '').trim().toLowerCase();
                            const download = (a.getAttribute('download') || '').toLowerCase();
                            const href = a.getAttribute('href') || '';

                            return text.includes('download snapshot') ||
                                   download.endsWith('.svg') ||
                                   href.startsWith('data:image/svg+xml');
                        });

                        if (!downloadLink) {
                            return 'download_link_not_found';
                        }

                        downloadLink.click();
                        return 'download_link_clicked';
                    }
                    """
                )

            download = download_info.value
            suggested_filename = download.suggested_filename

            temp_svg_path = out_path.with_suffix(".tmp.svg")
            download.save_as(str(temp_svg_path))
            svg_text = temp_svg_path.read_text(encoding="utf-8")
            temp_svg_path.unlink(missing_ok=True)

        except PlaywrightTimeoutError:
            clicked = "download_timeout_fallback_extract_data_uri"

            data_uri = popup.evaluate(
                """
                () => {
                    const links = Array.from(document.querySelectorAll('a'));
                    const downloadLink = links.find(a => {
                        const text = (a.innerText || a.textContent || '').trim().toLowerCase();
                        const download = (a.getAttribute('download') || '').toLowerCase();
                        const href = a.getAttribute('href') || '';

                        return text.includes('download snapshot') ||
                               download.endsWith('.svg') ||
                               href.startsWith('data:image/svg+xml');
                    });

                    if (!downloadLink) {
                        return null;
                    }

                    return downloadLink.getAttribute('href') || null;
                }
                """
            )

            svg_text = extract_svg_from_data_uri(data_uri) if data_uri else None

        browser.close()

    if not svg_text:
        raise RuntimeError(
            "Could not download or extract SVG from the Krona snapshot popup. "
            f"Snapshot trigger result: {triggered}; download action: {clicked}"
        )

    svg_start = svg_text.find("<svg")
    if svg_start == -1:
        preview = svg_text[:500].replace("\n", "\\n")
        raise RuntimeError(
            "Output was captured, but no <svg> tag was found. "
            f"First 500 characters were: {preview}"
        )

    svg_text = svg_text[svg_start:]

    if not args.no_white_background:
        svg_text = add_white_background(svg_text)

    svg_text = scale_svg_font_sizes(
        svg_text,
        scale=args.font_scale,
        fixed_px=args.font_size,
    )

    out_path.write_text(svg_text, encoding="utf-8")

    if png_path:
        svg_to_png_with_cairosvg(
            svg_path=out_path,
            png_path=png_path,
            output_width=args.png_width,
            output_height=args.png_height,
        )

    print(f"Triggered snapshot using: {triggered}")
    print(f"Download action: {clicked}")

    if suggested_filename:
        print(f"Downloaded suggested filename: {suggested_filename}")

    print(f"Wrote SVG: {out_path}")

    if png_path:
        print(f"Wrote PNG: {png_path}")


if __name__ == "__main__":
    main()
