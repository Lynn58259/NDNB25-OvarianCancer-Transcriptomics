"""Compare regenerated numerical tables with a separate reference checkout."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd


def compare(reference, generated, atol=1e-10, rtol=1e-8):
    rows = []
    files = sorted(p for p in (reference / 'results').rglob('*.csv')
                   if 'quality_control' not in p.relative_to(reference / 'results').parts)
    if not files:
        raise ValueError('No reference result CSV files found.')
    for ref in files:
        rel = ref.relative_to(reference)
        dest = generated / rel
        row = dict(file=rel.as_posix(), match=False, reason='', numeric_mismatches=0,
                   text_mismatches=0, max_absolute_difference=0., significance_changes=0)
        if not dest.is_file():
            row['reason'] = 'missing generated file'; rows.append(row); continue
        x = pd.read_csv(ref, keep_default_na=False)
        y = pd.read_csv(dest, keep_default_na=False)
        if x.shape != y.shape or list(x.columns) != list(y.columns):
            row['reason'] = 'shape or columns differ'; rows.append(row); continue
        for col in x:
            xn = pd.to_numeric(x[col], errors='coerce')
            yn = pd.to_numeric(y[col], errors='coerce')
            valid = xn.notna() & yn.notna()
            xx = xn[valid].to_numpy(dtype=float)
            yy = yn[valid].to_numpy(dtype=float)
            row['numeric_mismatches'] += int((~np.isclose(xx, yy, atol=atol, rtol=rtol)).sum())
            finite = np.isfinite(xx) & np.isfinite(yy)
            if finite.any():
                row['max_absolute_difference'] = max(row['max_absolute_difference'],
                                                     float(abs(xx[finite] - yy[finite]).max()))
            row['text_mismatches'] += int((x.loc[~valid, col].astype(str) !=
                                          y.loc[~valid, col].astype(str)).sum())
            if col in ('padj', 'FDR', 'BH_q_across_three_prespecified_tests'):
                row['significance_changes'] += int(((xn < .05) != (yn < .05)).sum())
        row['match'] = not (row['numeric_mismatches'] or row['text_mismatches'] or
                            row['significance_changes'])
        rows.append(row)
    return pd.DataFrame(rows)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--generated', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=Path('comparison.csv'))
    parser.add_argument('--atol', type=float, default=1e-10)
    parser.add_argument('--rtol', type=float, default=1e-8)
    args = parser.parse_args()
    if args.reference.resolve() == args.generated.resolve():
        parser.error('Reference and generated directories must be different.')
    result = compare(args.reference.resolve(), args.generated.resolve(), args.atol, args.rtol)
    result.to_csv(args.output, index=False)
    passed = int(result['match'].sum())
    print(f'{passed}/{len(result)} tables match; {int(result.significance_changes.sum())} significance changes.')
    raise SystemExit(0 if result['match'].all() else 1)
