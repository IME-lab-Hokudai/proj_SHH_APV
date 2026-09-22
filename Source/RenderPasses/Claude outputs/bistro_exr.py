import os
# OpenCV only reads/writes EXR if this is set BEFORE cv2 is imported.
os.environ["OPENCV_IO_ENABLE_OPENEXR"] = "1"

import cv2
import numpy as np
from PIL import Image, ImageDraw, ImageFont

try:
    import matplotlib.font_manager as fm
    DEFAULT_FONT_PATH = fm.findfont("DejaVu Sans")
except Exception:
    DEFAULT_FONT_PATH = None


# ============================================================
# USER SETTINGS
# ============================================================

# Linear HDR renders (EXR, all the same resolution), captured BEFORE tone mapping,
# with the measurement shader settings (no SH window, no clamp).
u64_render_path = "U64.exr"
hessian_render_path = "Hessian.exr"
egc_render_path = "EGC.exr"

# Optional mask (EXR or PNG, same resolution): nonzero = pixels to evaluate,
# e.g. the dynamic-object mask of the demo pass. None = all pixels with
# reference luminance > epsilon (same as before).
mask_path = None

output_image_path = "figure_zoom_compare_output.png"

resize_factor = 0.40
epsilon = 1e-5
max_heatmap_error = 10.0  # percent, top of the colorbar

# ---- Display transform (ONLY for the shown images, never for the metrics) ----
# exposure -> Reinhard x/(1+x) -> sRGB encode -> 8 bit. Keep identical for all
# panels and state it in the caption.
display_exposure = 1.0

# ---- PSNR ----
# Linear HDR has no natural peak value, so PSNR is computed on the tone-mapped
# display images (same transform as above, values in [0,1]).
# MAPE and RMSE are computed on linear luminance.

# Zoom box, defined by its center and size as fractions of the resized
# panel (new_w, new_h). (0,0) is top-left, (1,1) is bottom-right.
crop_center_frac = (0.28, 0.8)   # (x, y) center of the box
crop_size_frac = (0.15, 0.17)     # (width, height) of the box

# The zoomed crop is pasted directly on top of its own region, enlarged
# around the same center point, rather than off in a corner.
zoom_multiplier = 1.8        # how much bigger the overlay is than the crop box
zoom_border_color = (13, 79, 133)   # BGR
zoom_border_thickness = 2

# Layout
outer_margin = 30
panel_gap = 24
row_gap_extra = 12   # extra vertical padding around the row-2 header
header_height = 80
bottom_margin = 34

header_font_size = 36
metric_font_size = 24
colorbar_font_size = 20
caption_font_size = 18

colorbar_width = 26
colorbar_label_width = 60
colorbar_height_frac = 0.55  # fraction of the reference column height


# ============================================================
# HELPERS
# ============================================================

def get_font(size):
    if DEFAULT_FONT_PATH is not None:
        return ImageFont.truetype(DEFAULT_FONT_PATH, size)
    return ImageFont.load_default()


def load_exr(path):
    """Load a linear float EXR as float32 BGR (alpha dropped). Values are kept
    as stored: unbounded and possibly negative."""
    img = cv2.imread(path, cv2.IMREAD_UNCHANGED | cv2.IMREAD_ANYDEPTH | cv2.IMREAD_ANYCOLOR)
    if img is None:
        raise FileNotFoundError(
            f"Could not load EXR: {path} "
            "(check the path and that this OpenCV build supports OpenEXR)")
    img = img.astype(np.float32)
    if img.ndim == 2:
        img = np.repeat(img[:, :, None], 3, axis=2)
    return img[:, :, :3]


def load_mask(path, shape):
    m = cv2.imread(path, cv2.IMREAD_UNCHANGED | cv2.IMREAD_ANYDEPTH | cv2.IMREAD_ANYCOLOR)
    if m is None:
        raise FileNotFoundError(f"Could not load mask: {path}")
    if m.ndim == 3:
        m = m[:, :, 2] if m.shape[2] >= 3 else m[:, :, 0]  # red channel (BGR order)
    if m.shape[:2] != shape[:2]:
        raise ValueError("Mask resolution does not match the renders.")
    return m.astype(np.float32) > 0.5 * float(m.max() if m.max() > 0 else 1.0)


def to_display_bgr(img_linear):
    """Linear HDR -> 8-bit display: exposure, Reinhard, sRGB. Display only."""
    x = np.maximum(img_linear * display_exposure, 0.0)   # negative light can't be displayed
    x = x / (1.0 + x)                                      # Reinhard
    srgb = np.where(x <= 0.0031308, 12.92 * x, 1.055 * np.power(x, 1.0 / 2.4) - 0.055)
    return (np.clip(srgb, 0.0, 1.0) * 255.0 + 0.5).astype(np.uint8)


def to_display_float(img_linear):
    """Same display transform, as float in [0,1] (used for PSNR)."""
    return to_display_bgr(img_linear).astype(np.float32) / 255.0


def resize_display(img_linear, new_w, new_h):
    # Downsample in linear space (physically correct averaging), then map to display.
    small = cv2.resize(img_linear, (new_w, new_h), interpolation=cv2.INTER_AREA)
    return to_display_bgr(small)


def center_pad(img_uint8, target_w, target_h, bg=255):
    h, w = img_uint8.shape[:2]
    canvas = np.full((target_h, target_w, 3), bg, dtype=np.uint8)
    x_off = (target_w - w) // 2
    y_off = (target_h - h) // 2
    canvas[y_off:y_off + h, x_off:x_off + w] = img_uint8
    return canvas


def luminance_bgr(img):
    return (
        0.2126 * img[:, :, 2] +
        0.7152 * img[:, :, 1] +
        0.0722 * img[:, :, 0]
    )


def compute_metrics(test_img, ref_img, mask):
    """MAPE and RMSE on LINEAR luminance; PSNR on the tone-mapped display images."""
    test_lum = luminance_bgr(test_img)
    ref_lum = luminance_bgr(ref_img)

    difference = test_lum - ref_lum
    absolute_error = np.abs(difference)

    ape = (
        absolute_error /
        np.maximum(np.abs(ref_lum), epsilon)
    ) * 100.0

    mape = float(np.mean(ape[mask]))
    rmse = float(np.sqrt(np.mean(difference[mask] ** 2)))

    # PSNR in display space (tone-mapped luminance in [0,1]).
    test_disp = luminance_bgr(to_display_float(test_img))
    ref_disp = luminance_bgr(to_display_float(ref_img))
    mse_disp = float(np.mean((test_disp[mask] - ref_disp[mask]) ** 2))
    psnr = 100.0 if mse_disp <= 0.0 else float(10.0 * np.log10(1.0 / mse_disp))

    return {"ape": ape, "mape": mape, "rmse": rmse, "psnr": psnr}


def make_text_panel(width, height, text, font_size=22):
    img = Image.new("RGB", (width, height), (255, 255, 255))
    draw = ImageDraw.Draw(img)
    font = get_font(font_size)

    lines = text.split("\n")
    boxes = [draw.textbbox((0, 0), line, font=font) for line in lines]
    line_heights = [box[3] - box[1] for box in boxes]

    total_h = sum(line_heights) + 4 * (len(lines) - 1)
    y = (height - total_h) // 2

    for line, box, line_height in zip(lines, boxes, line_heights):
        text_width = box[2] - box[0]
        x = (width - text_width) // 2
        draw.text((x, y), line, font=font, fill=(0, 0, 0))
        y += line_height + 4

    return cv2.cvtColor(np.array(img), cv2.COLOR_RGB2BGR)


def draw_text_bgr(img_bgr, text, xy, font_size=22, color=(255, 255, 255)):
    img_rgb = cv2.cvtColor(img_bgr, cv2.COLOR_BGR2RGB)
    pil_img = Image.fromarray(img_rgb)
    draw = ImageDraw.Draw(pil_img)
    draw.text(xy, text, font=get_font(font_size), fill=color)
    return cv2.cvtColor(np.array(pil_img), cv2.COLOR_RGB2BGR)


def add_zoom_overlay(panel_bgr, crop_box):
    h, w = panel_bgr.shape[:2]
    x0, y0, x1, y1 = crop_box
    crop = panel_bgr[y0:y1, x0:x1]
    if crop.size == 0:
        return panel_bgr.copy()

    crop_w, crop_h = x1 - x0, y1 - y0
    cx, cy = (x0 + x1) / 2.0, (y0 + y1) / 2.0

    zoom_w = max(1, int(round(crop_w * zoom_multiplier)))
    zoom_h = max(1, int(round(crop_h * zoom_multiplier)))
    zoom_resized = cv2.resize(crop, (zoom_w, zoom_h), interpolation=cv2.INTER_NEAREST)

    zx0 = int(round(cx - zoom_w / 2.0))
    zy0 = int(round(cy - zoom_h / 2.0))
    zx0 = max(0, min(zx0, w - zoom_w))
    zy0 = max(0, min(zy0, h - zoom_h))
    zx1, zy1 = zx0 + zoom_w, zy0 + zoom_h

    out = panel_bgr.copy()
    out[zy0:zy1, zx0:zx1] = zoom_resized
    cv2.rectangle(out, (zx0, zy0), (zx1, zy1), zoom_border_color, zoom_border_thickness)

    return out


def make_error_heatmap(ape, mask, new_w, new_h, mape=None, rmse=None, psnr=None):
    normalized_error = np.clip(ape / max_heatmap_error, 0.0, 1.0)
    normalized_error = (normalized_error * 255.0).astype(np.uint8)

    heatmap = cv2.applyColorMap(normalized_error, cv2.COLORMAP_JET)
    heatmap[~mask] = [0, 0, 0]
    heatmap = cv2.resize(heatmap, (new_w, new_h), interpolation=cv2.INTER_NEAREST)

    metric_lines = []
    if mape is not None:
        metric_lines.append(f"MAPE  {mape:.3f}%")
    if rmse is not None:
        metric_lines.append(f"RMSE  {rmse:.4f}")
    if psnr is not None:
        metric_lines.append(f"PSNR  {psnr:.2f} dB (tm)")

    if metric_lines:
        img_rgb = cv2.cvtColor(heatmap, cv2.COLOR_BGR2RGB)
        pil_img = Image.fromarray(img_rgb)
        draw = ImageDraw.Draw(pil_img)
        font = get_font(metric_font_size)

        x, y = 14, 14
        padding_x, padding_y, line_gap = 10, 7, 5

        boxes = [draw.textbbox((0, 0), line, font=font) for line in metric_lines]
        text_width = max(box[2] - box[0] for box in boxes)
        heights = [box[3] - box[1] for box in boxes]
        text_height = sum(heights) + line_gap * (len(metric_lines) - 1)

        draw.rectangle(
            (x - padding_x, y - padding_y, x + text_width + padding_x, y + text_height + padding_y),
            fill=(0, 0, 0),
        )

        yy = y
        for line, line_height in zip(metric_lines, heights):
            draw.text((x, yy), line, font=font, fill=(255, 255, 255))
            yy += line_height + line_gap

        heatmap = cv2.cvtColor(np.array(pil_img), cv2.COLOR_RGB2BGR)

    return heatmap


def make_color_bar(width, height):
    values = np.linspace(255, 0, height, dtype=np.uint8).reshape(height, 1)
    values = np.repeat(values, width, axis=1)
    return cv2.applyColorMap(values, cv2.COLORMAP_JET)


def paste(canvas, img, x, y):
    h, w = img.shape[:2]
    canvas[y:y + h, x:x + w] = img


# ============================================================
# LOAD + COMPUTE
# ============================================================

u64_img = load_exr(u64_render_path)
hessian_img = load_exr(hessian_render_path)
egc_img = load_exr(egc_render_path)

if hessian_img.shape != u64_img.shape or egc_img.shape != u64_img.shape:
    raise ValueError("All three renders must have the same resolution.")

u64_lum = luminance_bgr(u64_img)
if mask_path is not None:
    mask = load_mask(mask_path, u64_img.shape) & (np.abs(u64_lum) > epsilon)
else:
    mask = np.abs(u64_lum) > epsilon

metrics_hessian = compute_metrics(hessian_img, u64_img, mask)
metrics_egc = compute_metrics(egc_img, u64_img, mask)

print(f"Hessian-only : MAPE {metrics_hessian['mape']:.3f}%  RMSE {metrics_hessian['rmse']:.5f}  "
      f"PSNR(tm) {metrics_hessian['psnr']:.2f} dB")
print(f"Hessian+EGC  : MAPE {metrics_egc['mape']:.3f}%  RMSE {metrics_egc['rmse']:.5f}  "
      f"PSNR(tm) {metrics_egc['psnr']:.2f} dB")

new_w = int(u64_img.shape[1] * resize_factor)
new_h = int(u64_img.shape[0] * resize_factor)

crop_w_px = int(new_w * crop_size_frac[0])
crop_h_px = int(new_h * crop_size_frac[1])
crop_cx_px = int(new_w * crop_center_frac[0])
crop_cy_px = int(new_h * crop_center_frac[1])

crop_x0 = max(0, min(crop_cx_px - crop_w_px // 2, new_w - crop_w_px))
crop_y0 = max(0, min(crop_cy_px - crop_h_px // 2, new_h - crop_h_px))
crop_box = (crop_x0, crop_y0, crop_x0 + crop_w_px, crop_y0 + crop_h_px)

u64_render_panel = resize_display(u64_img, new_w, new_h)
u64_render_panel = add_zoom_overlay(u64_render_panel, crop_box)

hessian_render_panel = resize_display(hessian_img, new_w, new_h)
egc_render_panel = resize_display(egc_img, new_w, new_h)

hessian_render_panel = add_zoom_overlay(hessian_render_panel, crop_box)
egc_render_panel = add_zoom_overlay(egc_render_panel, crop_box)

hessian_heatmap = make_error_heatmap(
    metrics_hessian["ape"], mask, new_w, new_h,
    mape=metrics_hessian["mape"], rmse=metrics_hessian["rmse"], psnr=metrics_hessian["psnr"],
)
egc_heatmap = make_error_heatmap(
    metrics_egc["ape"], mask, new_w, new_h,
    mape=metrics_egc["mape"], rmse=metrics_egc["rmse"], psnr=metrics_egc["psnr"],
)


# ============================================================
# LAYOUT
# ============================================================

col1_x = outer_margin
col2_x = col1_x + new_w + panel_gap
col3_x = col2_x + new_w + panel_gap
colorbar_x = col3_x + new_w + panel_gap

total_width = colorbar_x + colorbar_width + colorbar_label_width + outer_margin

y0 = outer_margin

# Row 1: image, then its label below it
row1_top = y0
row1_bottom = row1_top + new_h
row1_label_top = row1_bottom

# Row 2: same pattern, offset below row 1's label
row2_top = row1_label_top + header_height + row_gap_extra
row2_bottom = row2_top + new_h
row2_label_top = row2_bottom

# Reference column spans the same vertical extent as rows 1+2 combined,
# with its own label below it (aligned with row 2's label row)
ref_top = row1_top
ref_bottom = row2_bottom
ref_height = ref_bottom - ref_top
ref_label_top = ref_bottom

caption_top = row2_label_top + header_height + 20
total_height = caption_top + caption_font_size + bottom_margin

canvas = np.full((total_height, total_width, 3), 255, dtype=np.uint8)

# Reference column (tall letterboxed render, label below)
paste(canvas, center_pad(u64_render_panel, new_w, ref_height), col1_x, ref_top)
paste(canvas, make_text_panel(new_w, header_height, "U64 reference render", header_font_size), col1_x, ref_label_top)

# Row 1: Hessian-only
paste(canvas, hessian_render_panel, col2_x, row1_top)
paste(canvas, hessian_heatmap, col3_x, row1_top)
paste(canvas, make_text_panel(new_w, header_height, "Hessian-only render", header_font_size), col2_x, row1_label_top)
paste(canvas, make_text_panel(new_w, header_height, "Hessian-only error map", header_font_size), col3_x, row1_label_top)

# Row 2: Hessian + EGC
paste(canvas, egc_render_panel, col2_x, row2_top)
paste(canvas, egc_heatmap, col3_x, row2_top)
paste(canvas, make_text_panel(new_w, header_height, "Hessian+EGC render", header_font_size), col2_x, row2_label_top)
paste(canvas, make_text_panel(new_w, header_height, "Hessian+EGC error map", header_font_size), col3_x, row2_label_top)

# Shared colorbar, shorter than and centered within the reference column
colorbar_height = int(ref_height * colorbar_height_frac)
colorbar_top = ref_top + (ref_height - colorbar_height) // 2
colorbar_bottom = colorbar_top + colorbar_height

color_bar = make_color_bar(colorbar_width, colorbar_height)
paste(canvas, color_bar, colorbar_x, colorbar_top)

canvas = draw_text_bgr(
    canvas, f"{max_heatmap_error:.0f}%",
    (colorbar_x + colorbar_width + 8, colorbar_top - 4),
    font_size=colorbar_font_size, color=(0, 0, 0),
)
canvas = draw_text_bgr(
    canvas, "Error",
    (colorbar_x + colorbar_width + 8, colorbar_top + colorbar_height // 2 - 12),
    font_size=colorbar_font_size, color=(0, 0, 0),
)
canvas = draw_text_bgr(
    canvas, "0%",
    (colorbar_x + colorbar_width + 8, colorbar_bottom - 24),
    font_size=colorbar_font_size, color=(0, 0, 0),
)

# Caption: document the display transform and the metric spaces.
canvas = draw_text_bgr(
    canvas,
    f"Renders: exposure {display_exposure:g}, Reinhard, sRGB (display only). "
    "MAPE/RMSE on linear luminance; PSNR on tone-mapped luminance.",
    (col1_x, caption_top),
    font_size=caption_font_size, color=(80, 80, 80),
)

cv2.imwrite(output_image_path, canvas)
print(f"Saved figure: {output_image_path}")
