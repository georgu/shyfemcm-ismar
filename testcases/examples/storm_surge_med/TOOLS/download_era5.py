#!/usr/bin/env python3
"""
Scarica il forcing atmosferico ERA5 (vento 10 m + pressione al livello del mare)
per il dominio mediterraneo, ottobre 2018 + un margine in coda.

Strategia: due richieste separate (ottobre intero, novembre 1-2) unite lungo
l'asse temporale, cosi' NON si scarica tutto novembre.

Prerequisiti (una volta sola):
  1) account sul Climate Data Store (CDS)
  2)  ~/.cdsapirc  con url + key (pagina "How to use the API" del CDS)
  3)  pip install cdsapi
  4)  cdo installato per il merge (sudo apt install cdo). In alternativa vedi
      la nota in fondo per unire con xarray.

Uso:
  python download_era5.py
"""

import os
import subprocess
import sys
import cdsapi

# ----------------------------------------------------------------------------
# CONFIGURAZIONE  (cambia qui per altri periodi / domini)
# ----------------------------------------------------------------------------
DATASET   = "reanalysis-era5-single-levels"
VARIABLES = [
    "10m_u_component_of_wind",
    "10m_v_component_of_wind",
    "mean_sea_level_pressure",
]
# Area del dominio + ~0.5 gradi di margine.  Ordine ERA5: [Nord, Ovest, Sud, Est]
AREA = [48, -11, 28, 41]
HOURS = [f"{h:02d}:00" for h in range(24)]   # tutte le 24 ore

# Blocco 1: ottobre 2018, giorni 1-31
REQ_OCT = {
    "product_type": "reanalysis",
    "variable": VARIABLES,
    "year": "2018",
    "month": "10",
    "day": [f"{d:02d}" for d in range(1, 32)],
    "time": HOURS,
    "area": AREA,
    "data_format": "netcdf",   # sul CDS "legacy" la chiave era "format"
}

# Blocco 2: novembre 2018, solo giorni 1-2 (il margine per itend = 1 nov 23:00)
REQ_NOV = {
    "product_type": "reanalysis",
    "variable": VARIABLES,
    "year": "2018",
    "month": "11",
    "day": ["01", "02"],
    "time": HOURS,
    "area": AREA,
    "data_format": "netcdf",
}

FILE_OCT = "era5_oct2018.nc"
FILE_NOV = "era5_nov2018_12.nc"
FILE_OUT = "era5_wind_201810.nc"

# ----------------------------------------------------------------------------
# DOWNLOAD
# ----------------------------------------------------------------------------
client = cdsapi.Client()


def scarica(request, filename, etichetta):
    """Scarica solo se il file non c'e' gia' (e non e' vuoto)."""
    if os.path.isfile(filename) and os.path.getsize(filename) > 0:
        print(f">> {filename} gia' presente ({os.path.getsize(filename)/1e6:.1f} MB)"
              f" -- salto il download di {etichetta}.")
        return
    print(f">> Scarico {etichetta} ...")
    client.retrieve(DATASET, request, filename)


scarica(REQ_OCT, FILE_OCT, "ottobre 2018")
scarica(REQ_NOV, FILE_NOV, "1-2 novembre 2018")

# ----------------------------------------------------------------------------
# MERGE lungo il tempo
# cdo mergetime e' robusto: gestisce da solo il nome della coordinata temporale
# (time / valid_time) che il nuovo CDS a volte cambia.
# ----------------------------------------------------------------------------
if os.path.isfile(FILE_OUT) and os.path.getsize(FILE_OUT) > 0:
    print(f">> {FILE_OUT} gia' presente -- niente da fare.")
    sys.exit(0)

print(f">> Unisco i due file in {FILE_OUT} ...")
try:
    subprocess.run(["cdo", "mergetime", FILE_OCT, FILE_NOV, FILE_OUT], check=True)
    print(f">> Fatto: {FILE_OUT}")
except (FileNotFoundError, subprocess.CalledProcessError):
    print("!! cdo non disponibile o merge fallito.")
    print(f"!! I due file grezzi restano qui: {FILE_OCT} e {FILE_NOV}")
    print("!! Uniscili a mano con:  cdo mergetime "
          f"{FILE_OCT} {FILE_NOV} {FILE_OUT}")
    sys.exit(1)
