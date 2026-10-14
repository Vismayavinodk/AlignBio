from .models import AlignmentResult


def format_alignment(result: AlignmentResult, line_width: int = 60) -> str:
    """Formats an alignment for terminal or text display with match/mismatch indicators."""
    s1, s2 = result.seq1_aligned, result.seq2_aligned
    lines = []

    # Header summary
    lines.append(f"Alignment Type: {result.alignment_type}")
    lines.append(f"Matrix: {result.matrix_name} | Gap Open: {result.gap_open} | Gap Extend: {result.gap_extend}")
    lines.append(f"Score: {result.score:.1f} | Length: {result.length} | Identity: {result.identity:.1f}%")
    lines.append("-" * min(line_width + 15, 80))

    for chunk_start in range(0, len(s1), line_width):
        chunk1 = s1[chunk_start : chunk_start + line_width]
        chunk2 = s2[chunk_start : chunk_start + line_width]

        markers = []
        for c1, c2 in zip(chunk1, chunk2):
            if c1 == "-" or c2 == "-":
                markers.append(" ")
            elif c1.upper() == c2.upper():
                markers.append("|")
            else:
                markers.append(":")

        lines.append(f"Seq1  {chunk_start+1:5d}  {chunk1}")
        lines.append(f"             {''.join(markers)}")
        lines.append(f"Seq2  {chunk_start+1:5d}  {chunk2}\n")

    return "\n".join(lines)


def alignment_to_html(result: AlignmentResult) -> str:
    """Generates an HTML snippet visualizing the alignment with CSS highlighting."""
    s1, s2 = result.seq1_aligned, result.seq2_aligned

    cells_top = []
    cells_mid = []
    cells_bot = []

    for c1, c2 in zip(s1, s2):
        if c1 == "-" or c2 == "-":
            bg = "#f3f4f6"
            marker = "&nbsp;"
        elif c1.upper() == c2.upper():
            bg = "#dcfce7"
            marker = "|"
        else:
            bg = "#fef9c3"
            marker = ":"

        cells_top.append(f'<span style="background:{bg}; font-family:monospace; padding:1px 3px;">{c1}</span>')
        cells_mid.append(f'<span style="font-family:monospace; padding:1px 3px; color:#6b7280;">{marker}</span>')
        cells_bot.append(f'<span style="background:{bg}; font-family:monospace; padding:1px 3px;">{c2}</span>')

    html = f"""<div style="font-family: sans-serif; border: 1px solid #e5e7eb; border-radius: 8px; padding: 16px; background: #ffffff;">
  <h4 style="margin: 0 0 8px 0;">{result.alignment_type.title()} Alignment (Score: {result.score:.1f}, Identity: {result.identity:.1f}%)</h4>
  <div style="overflow-x: auto; white-space: nowrap; line-height: 1.4;">
    <div>{''.join(cells_top)}</div>
    <div>{''.join(cells_mid)}</div>
    <div>{''.join(cells_bot)}</div>
  </div>
</div>"""
    return html
