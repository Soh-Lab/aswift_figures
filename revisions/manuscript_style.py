"""Reuse the original notebooks' definitions without running their data-loading cells."""
from .common import ROOT, METHODS, LABELS, COLORS
import ast
import json
import matplotlib.pyplot as plt


def definitions(name, stop_assignment):
    path = ROOT / 'figures' / f'{name}.ipynb'
    notebook = json.loads(path.read_text())
    source = '\n\n'.join(''.join(c['source']) for c in notebook['cells'] if c['cell_type'] == 'code')
    tree = ast.parse(source)
    body = []
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == stop_assignment for t in node.targets):
            break
        body.append(node)
    else:
        raise ValueError(f'Missing definition boundary {stop_assignment} in {path}')
    # Prevent style settings from an earlier figure (e.g. missing spines) leaking in.
    plt.rcdefaults()
    namespace = {'__name__': f'revisions.style_{name}'}
    exec(compile(ast.Module(body=body, type_ignores=[]), str(path), 'exec'), namespace)
    return namespace


def figure2(methods=None, labels=None):
    ns = definitions('FIG_2', 'panel_summary_path')
    if methods is not None:
        template = ns['METHOD_STYLE']['aswift'].copy()
        # Mutate in place: original helper defaults reference these list/dict objects.
        ns['METHOD_ORDER'][:] = methods
        ns['METHOD_LABELS'].clear()
        ns['METHOD_LABELS'].update(dict(zip(methods, labels or methods)))
        ns['METHOD_STYLE'].clear()
        palette = COLORS + ['#9467bd', '#8c564b']
        for i, method in enumerate(methods):
            ns['METHOD_STYLE'][method] = dict(template, color=palette[i % len(palette)])
        for key, field in [('METHOD_COLORS', 'color'), ('METHOD_MARKERS', 'marker'), ('METHOD_LINESTYLES', 'linestyle')]:
            ns[key].clear()
            ns[key].update({m: ns['METHOD_STYLE'][m][field] for m in methods})
    ns['set_figure_style']()
    return ns


def figure4():
    return definitions('FIG_4', 'methods')
