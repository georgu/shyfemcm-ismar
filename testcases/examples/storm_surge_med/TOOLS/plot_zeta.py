#!/usr/bin/env python3
"""
plot_zeta.py -- water level (zeta) vs time ai 4 mareografi, con taglio dello spin-up.

Uso (da dentro sims/all o sims/tide):
    python3 ../../plot_zeta.py
    python3 ../../plot_zeta.py --skip 7

Uso (da sims/, senza spostarsi):
    python3 ../plot_zeta.py --dir all
    python3 ../plot_zeta.py --dir tide

L'etichetta del run e' dedotta dal nome della cartella; sovrascrivibile con --label.
Il cold start parte da livello zero: i primi giorni sono transitorio numerico,
non fisica. 5 giorni bastano; Vaia (29 ott) resta abbondantemente dentro.
"""

import os
import glob
import argparse
from datetime import datetime, timedelta
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

STATIONS = {'1': 'Trieste', '2': 'Venezia', '3': 'Ravenna', '4': 'Ancona'}
ORDER = ['Trieste', 'Venezia', 'Ravenna', 'Ancona']


def load(path):
    t, z = [], []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            p = line.split()
            if len(p) < 2:
                continue
            try:
                tt = datetime.strptime(p[0], '%Y-%m-%d::%H:%M:%S')
                zz = float(p[1])
            except ValueError:
                continue
            t.append(tt)
            z.append(zz)
    return t, z


def trim(t, z, cutoff):
    tt = [a for a in t if a >= cutoff]
    zz = [b for a, b in zip(t, z) if a >= cutoff]
    return tt, zz


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dir', default='.', help='folder with files zeta.2d.* (default: cwd)')
    ap.add_argument('--skip', type=float, default=5.0,
                    help='spin-up days to be excluded (default 5)')
    ap.add_argument('--label', default=None,
                    help='label of run (default: folder name)')
    args = ap.parse_args()

    if not os.path.isdir(args.dir):
        print('Folder not found: %s' % args.dir)
        return

    files = sorted(glob.glob(os.path.join(args.dir, 'zeta.2d.*')))
    if not files:
        print('No zeta.2d.* file in %s' % os.path.abspath(args.dir))
        return

    label = args.label or os.path.basename(os.path.abspath(args.dir))

    all_first = None
    for path in files:
        t, _ = load(path)
        if t:
            all_first = t[0] if all_first is None else min(all_first, t[0])
    if all_first is None:
        print('Files are empty or cannot be read.')
        return
    cutoff = all_first + timedelta(days=args.skip)

    print('run         : %s' % label)
    print('data start : %s' % all_first.strftime('%Y-%m-%d %H:%M'))
    print('spin-up cut-off: %g days  ->  plot from %s'
          % (args.skip, cutoff.strftime('%Y-%m-%d %H:%M')))
    print()
    print('%-10s %8s %8s %10s   %-20s' %
          ('station', 'min[m]', 'max[m]', 'excurs[m]', 'hour of max'))
    print('   (stats about the chosen time window, excluding spin-up)')
    print('-' * 66)

    data = {}
    for path in files:
        idx = os.path.basename(path).split('.')[-1]
        name = STATIONS.get(idx, os.path.basename(path))
        t, z = load(path)
        t, z = trim(t, z, cutoff)
        if not z:
            print('%-10s  (null after cut)' % name)
            continue
        zmin, zmax = min(z), max(z)
        tmax = t[z.index(zmax)]
        exc = zmax - zmin
        flag = '  <-- suspect excursion (<5 cm)' if exc < 0.05 else ''
        print('%-10s %8.3f %8.3f %10.3f   %s%s' %
              (name, zmin, zmax, exc, tmax.strftime('%Y-%m-%d %H:%M'), flag))
        data[name] = (t, z)
    print()

    if not data:
        return

    fig, ax = plt.subplots(figsize=(13, 6))
    for name in ORDER:
        if name in data:
            t, z = data[name]
            ax.plot(t, z, lw=0.8, label=name)

    vaia = datetime(2018, 10, 29)
    if vaia >= cutoff:
        ax.axvline(vaia, color='0.6', ls='--', lw=0.8, zorder=0)
        ax.text(vaia, ax.get_ylim()[1], ' Vaia (29 Oct)',
                fontsize=8, color='0.4', va='top')

    ax.set_xlabel('time')
    ax.set_ylabel('zeta [m] (in respect to reference model sea level)')
    ax.set_title('Sea level at Adriatic gauges - run "%s" '
                 '(excluding %g days of spin-up)' % (label, args.skip))
    ax.legend(loc='upper left', fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%d %b'))
    ax.xaxis.set_major_locator(mdates.DayLocator(interval=2))
    fig.autofmt_xdate()

    out = os.path.join(args.dir, 'zeta_%s.png' % label)
    fig.savefig(out, dpi=150, bbox_inches='tight')
    print('Figure saved in: %s' % out)


if __name__ == '__main__':
    main()
