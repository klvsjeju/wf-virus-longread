"""Local report header theme based on the user's supplied reference artwork."""
import base64
from pathlib import Path
from dominate import tags as html
from dominate.util import raw

LOGO_LINK = 'https://github.com/klvsjeju/wf-virus-longread'

def descendants(node):
    yield node
    for child in getattr(node, 'children', []):
        if hasattr(child, 'children'):
            yield from descendants(child)

def apply_branding(report):
    """Style only the header; retain navigation controls and dynamic metadata."""
    artwork = Path(__file__).with_name('report_assets') / 'header-reference.png'
    uri = 'data:image/png;base64,' + base64.b64encode(artwork.read_bytes()).decode('ascii')
    # Viewport onto the circular logo in the original, unmodified reference PNG.
    # Keeping the asset embedded makes copied/downloaded HTML self-contained.
    for node in descendants(report):
        attrs = getattr(node, 'attributes', {})
        if attrs.get('href') == 'https://labs.epi2me.io/':
            node['href'] = LOGO_LINK
            node['aria-label'] = 'wf-virus-longread GitHub repository'
            node['class'] = 'workflow-logo-link me-auto text-decoration-none'
            node.children.clear()
            node.add(html.span(role='img', aria_label='Circular bird logo', cls='workflow-logo'))
            break
    else:
        raise ValueError('Report header logo link was not found; check LabsReport layout')
    badges = [node for node in descendants(report.banner)
              if 'badge' in getattr(node, 'attributes', {}).get('class', '').split()]
    if badges:
        badges[0].children.clear()
        badges[0].add('Research only')
    css = '''
    :root { --workflow-red: #C00000; --workflow-badge: #B00020; --workflow-ink: #042433; }
    header { background: #fff; border: 7px solid var(--workflow-red); }
    header nav.fixed-top.bg-dark { background-color: #fff !important; border: 7px solid var(--workflow-red); border-bottom: 0; }
    header nav .container { padding: 12px 24px !important; }
    .workflow-logo-link { display: inline-flex; align-items: center; }
    .workflow-logo { display: block; width: 76px; height: 76px; background-image: url("ARTWORK");
      background-repeat: no-repeat; background-size: 1534.06px 469.68px; background-position: -152px -91.2px; }
    .labs-banner.bg-dark { background-color: #fff !important; margin-bottom: 0 !important;
      padding: 34px 24px 48px !important; }
    .labs-banner .report-title { color: #000 !important; border: 0 !important; }
    .labs-banner h1 { color: #000; font-weight: 750; font-size: clamp(1.5rem, 3vw, 4.5rem); line-height: 1.15; margin: 12px 0 32px; }
    .labs-banner .report-title > p { display: none; }
    .labs-banner .badge { color: #fff !important; background-color: var(--workflow-badge) !important;
      border: 2px solid var(--workflow-ink); border-radius: 9px !important;
      font-size: 15px !important; font-weight: 400; padding: 12px 22px !important; }
    .labs-banner .badge:nth-child(2) { background-color: #CD1C18 !important; }
    .labs-banner .badge:nth-child(3) { background-color: #DC1C13 !important; }
    header nav #dropdownMenuButton { background-color: var(--workflow-badge); color: #fff;
      border: 2px solid var(--workflow-ink); border-radius: 9px; padding: 10px 18px; }
    header nav #dropdownMenuButton:hover, header nav #dropdownMenuButton:focus { background-color: #95001d; }
    @media (max-width: 600px) {
      header nav .container { padding: 6px 12px !important; }
      .labs-banner.bg-dark { padding: 24px 16px 30px !important; }
      .labs-banner .badge { padding: 8px 12px !important; font-size: 13px !important; }
    }
    '''.replace('ARTWORK', uri)
    report.banner.add(html.style(raw(css)))
