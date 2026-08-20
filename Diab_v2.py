#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 00:38:24 2026

@author: julienarmand
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LibreLinkUp — Glycemie + Pente WLS + Notifications ntfy + Watchdog healthchecks.io
Toutes les valeurs sont configurees directement ci-dessous.
"""

import os
import time
import hashlib
import traceback
from pathlib import Path

import requests
import numpy as np
import pandas as pd
from dotenv import load_dotenv

# ══════════════════════════════════════════════════════════════
# SECRETS — read from a .env file sitting next to this script.
# That file is git-ignored, so real values never reach the repository.
# Copy .env.example to .env and fill it in to run this script.
# ══════════════════════════════════════════════════════════════
load_dotenv(Path(__file__).resolve().parent / ".env")

# ── LibreLinkUp ───────────────────────────────────────────────
EMAIL    = os.environ.get("LIBRELINKUP_EMAIL", "")
PASSWORD = os.environ.get("LIBRELINKUP_PASSWORD", "")
COUNTRY  = os.environ.get("LIBRELINKUP_COUNTRY", "CA")

# ── ntfy and healthchecks.io ──────────────────────────────────
# The ntfy topic is a shared secret: anyone who knows it can read the
# glucose notifications. Same for the healthchecks.io UUID.
NTFY_TOPIC = os.environ.get("NTFY_TOPIC", "")
NTFY_URL   = os.environ.get("NTFY_URL", "https://ntfy.sh")
HC_UUID    = os.environ.get("HEALTHCHECKS_UUID", "")
HC_URL     = f"https://hc-ping.com/{HC_UUID}"

# Fail early with a readable message rather than a login error further down.
_missing = [name for name, value in (("LIBRELINKUP_EMAIL",    EMAIL),
                                     ("LIBRELINKUP_PASSWORD", PASSWORD),
                                     ("NTFY_TOPIC",           NTFY_TOPIC)) if not value]
if _missing:
    raise SystemExit(
        "Missing in .env: " + ", ".join(_missing) +
        "\nCopy .env.example to .env and fill in your own values."
    )

# The watchdog is optional: without a UUID it is simply switched off.
HC_ACTIF = bool(HC_UUID)
if not HC_ACTIF:
    print("[i] HEALTHCHECKS_UUID not set — watchdog disabled.")

# ── Nom de CETTE machine ──────────────────────────────────────
MACHINE = os.environ.get("MACHINE_NAME", "PC-A")

# ── Seuils slops (mmol/L/h) ───────────────────────────────────
SEUIL_HAUT =  6.0
SEUIL_BAS  = -4.5

# ── Seuils glycemie (mmol/L) ──────────────────────────────────
HYPO_SEUIL  = 5.0
HYPER_SEUIL = 14.0

# ── Intervals (seconds) ────────────────────────────────────
INTERVALLE     = 300   # entre chaque lecture (5 min)
SILENCE_ALERTE = 900   # entre deux alertes du meme type (15 min)

# ── Pente WLS ─────────────────────────────────────────────────
WLS_WINDOW = 8
WLS_DECAY  = 0.7


# ── healthchecks.io ────────────────────────────────────────────────
def hc_ping(suffixe="", payload=None):
    """
    Signale a healthchecks.io que le script tourne.
    suffixe : ""       -> tout va bien
              "/start" -> debut de cycle
              "/fail"  -> erreur (declenche l'alerte immediatement)
    """
    if not HC_ACTIF or HC_UUID.startswith("AAAAAAAA"):
        return
    try:
        requests.post(HC_URL + suffixe,
                      data=(payload or "").encode("utf-8"),
                      timeout=10)
    except Exception as e:
        print(f"  [hc FAIL] {e}")


# ── ntfy ──────────────────────────────────────────────────────────
def notif(titre, message, priorite="default", tags=""):
    """
    Envoie une notification ntfy.
    priorite : min | low | default | high | urgent
    """
    try:
        r = requests.post(
            f"{NTFY_URL}/{NTFY_TOPIC}",
            data=message.encode("utf-8"),
            headers={
                "Title":    f"[{MACHINE}] {titre}",
                "Priority": priorite,
                "Tags":     tags,
            },
            timeout=10,
        )
        ok = r.status_code == 200
        print(f"  [ntfy {'OK' if ok else 'ERR ' + str(r.status_code)}] {titre}")
    except Exception as e:
        print(f"  [ntfy FAIL] {e}")


# ── Pente WLS  ──────────────────────────────────────────────
def pente_wls(glucose_series, window=WLS_WINDOW, decay=WLS_DECAY):
    """Regression lineaire ponderee causale — aucun point futur utilise."""
    vals   = glucose_series.values.astype(float)
    result = np.full(len(vals), np.nan)
    w   = np.array([decay ** (window - 1 - i) for i in range(window)])
    w  /= w.sum()
    x   = np.arange(window, dtype=float)
    wx  = float(np.dot(w, x))
    wxx = float(np.dot(w, x ** 2))
    for i in range(window - 1, len(vals)):
        y = vals[i - window + 1 : i + 1]
        if np.any(np.isnan(y)):
            continue
        wy  = float(np.dot(w, y))
        wxy = float(np.dot(w, x * y))
        d   = wxx - wx ** 2
        if abs(d) > 1e-10:
            result[i] = (wxy - wx * wy) / d * 12
    return pd.Series(result, index=glucose_series.index)


def pente_ok(p):
    try:
        return p is not None and not np.isnan(float(p))
    except (TypeError, ValueError):
        return False


def fleche(p):
    if not pente_ok(p):
        return "?"
    p = float(p)
    if p >  3: return "↑↑"
    if p >  1: return "↑"
    if p < -3: return "↓↓"
    if p < -1: return "↓"
    return "→"


def to_mmol(val):
    val = float(val)
    return round(val / 18.018, 1) if val > 30 else val


# ── Client LibreLinkUp ─────────────────────────────────────────────
class LibreLinkUpClient:

    def __init__(self, email, password, country="CA"):
        self.email    = email
        self.password = password
        if country.upper() == "US":
            self.base_url = "https://api-us.libreview.io"
        elif country.upper() == "EU":
            self.base_url = "https://api-eu.libreview.io"
        else:
            self.base_url = "https://api.libreview.io"
        self.headers = {
            "version":         "4.16.0",
            "product":         "llu.android",
            "Accept-Encoding": "gzip",
            "Content-Type":    "application/json",
            "A-Component":     "Premium",
            "User-Agent":      "Mozilla/5.0 (Linux; Android 13; LLU) AppleWebKit/537.36",
        }
        self.token      = None
        self.patient_id = None

    def login(self):
        r = requests.post(f"{self.base_url}/llu/auth/login",
                          json={"email": self.email, "password": self.password},
                          headers=self.headers)
        if r.status_code == 200:
            data = r.json()
            if data.get("status") == 0:
                self.token   = data["data"]["authTicket"]["token"]
                user_id      = data["data"]["user"]["id"]
                account_hash = hashlib.sha256(user_id.encode()).hexdigest()
                self.headers["Authorization"] = f"Bearer {self.token}"
                self.headers["Account-Id"]    = account_hash
                print("[OK] Connexion reussie.")
                return True
        print(f"[X] Echec connexion : {r.text}")
        return False

    def fetch_patient(self):
        r    = requests.get(f"{self.base_url}/llu/connections", headers=self.headers)
        data = r.json().get("data", [])
        if data:
            self.patient_id = data[0]["patientId"]
            nom = f"{data[0]['firstName']} {data[0]['lastName']}"
            print(f"[OK] Patient : {nom}")
            return True
        print(f"[!] Aucun patient : {r.text}")
        return False

    def get_graph_data(self):
        if not self.patient_id:
            return None, []
        r = requests.get(
            f"{self.base_url}/llu/connections/{self.patient_id}/graph",
            headers=self.headers)
        if r.status_code != 200:
            print(f"[X] Erreur API : {r.status_code}")
            return None, []
        d = r.json().get("data", {})
        return (d.get("connection", {}).get("glucoseMeasurement"),
                d.get("graphData", []))


# ── Boucle principale ──────────────────────────────────────────────
def main():
    print("=" * 58)
    print(f"  MONITEUR GLYCEMIE — {MACHINE}")
    print(f"  Topic ntfy   : {NTFY_URL}/{NTFY_TOPIC}")
    print(f"  Watchdog     : {'actif' if HC_ACTIF else 'desactive'}")
    print(f"  Pente        : alerte si <= {SEUIL_BAS} ou >= +{SEUIL_HAUT} mmol/L/h")
    print(f"  Glycemie     : hypo < {HYPO_SEUIL}  |  hyper > {HYPER_SEUIL}")
    print(f"  Intervalle   : {INTERVALLE//60} min")
    print("=" * 58 + "\n")

    if HC_ACTIF and HC_UUID.startswith("AAAAAAAA"):
        print("  [!] HC_UUID n'a pas ete remplace — watchdog inactif.\n")

    client = LibreLinkUpClient(EMAIL, PASSWORD, COUNTRY)
    if not client.login() or not client.fetch_patient():
        hc_ping("/fail", "Echec de connexion LibreLinkUp au demarrage")
        return

    derniere_alerte = {}

    def peut_alerter(t):
        return (time.time() - derniere_alerte.get(t, 0)) >= SILENCE_ALERTE

    def marquer(t):
        derniere_alerte[t] = time.time()

    notif("Moniteur demarre",
          f"Surveillance glycemie active\nTopic: {NTFY_TOPIC}",
          priorite="low", tags="white_check_mark")
    hc_ping("/start", f"Demarrage du moniteur sur {MACHINE}")

    print("  Surveillance en cours... (Ctrl+C pour arreter)\n")

    while True:
        try:
            actuelle, historique = client.get_graph_data()

            if actuelle is None:
                print("  [!] Pas de donnees — attente 1 min.")
                hc_ping("/fail", "Aucune donnee retournee par l'API LibreLinkUp")
                time.sleep(60)
                continue

            # Construire le DataFrame
            lignes = []
            for pt in historique:
                val = pt.get("Value") or pt.get("value")
                ts  = pt.get("Timestamp") or pt.get("timestamp")
                if val is not None and ts:
                    lignes.append({
                        "timestamp": pd.to_datetime(ts),
                        "glucose":   to_mmol(val),
                    })
            lignes.append({
                "timestamp": pd.to_datetime(actuelle["Timestamp"]),
                "glucose":   to_mmol(actuelle["Value"]),
            })

            df = (pd.DataFrame(lignes)
                  .drop_duplicates("timestamp")
                  .sort_values("timestamp")
                  .reset_index(drop=True))

            # Pente WLS causale
            df["pente"] = pente_wls(df["glucose"])
            dernier     = df.iloc[-1]
            glucose     = float(dernier["glucose"])
            pente       = float(dernier["pente"]) if pd.notna(dernier["pente"]) else None
            heure       = dernier["timestamp"].strftime("%H:%M")
            zone        = ("HYPO"  if glucose < HYPO_SEUIL  else
                           "HYPER" if glucose > HYPER_SEUIL else "OK")
            pente_str   = f"{pente:+.2f}" if pente_ok(pente) else "  ?"

            ligne_statut = (f"[{heure}]  {glucose:.1f} mmol/L  {fleche(pente)}  "
                            f"(pente {pente_str} mmol/L/h)  [{zone}]")
            print(f"  {ligne_statut}")

            # ── Cycle reussi : on rassure healthchecks.io ──
            hc_ping("", f"{MACHINE} — {ligne_statut}")

            # Alertes
            if pente_ok(pente) and pente <= SEUIL_BAS and peut_alerter("pente_bas"):
                notif("Descente rapide",
                      f"Pente : {pente:+.2f} mmol/L/h\nGlycemie : {glucose:.1f} mmol/L ({heure})",
                      priorite="high", tags="chart_with_downwards_trend")
                marquer("pente_bas")

            if pente_ok(pente) and pente >= SEUIL_HAUT and peut_alerter("pente_haut"):
                notif("Montee rapide",
                      f"Pente : {pente:+.2f} mmol/L/h\nGlycemie : {glucose:.1f} mmol/L ({heure})",
                      priorite="high", tags="chart_with_upwards_trend")
                marquer("pente_haut")

            if glucose < HYPO_SEUIL and peut_alerter("hypo"):
                notif("HYPOGLYCEMIE",
                      f"Glycemie : {glucose:.1f} mmol/L\nPente : {pente_str} mmol/L/h ({heure})",
                      priorite="urgent", tags="rotating_light")
                marquer("hypo")

            if glucose > HYPER_SEUIL and peut_alerter("hyper"):
                notif("HYPERGLYCEMIE",
                      f"Glycemie : {glucose:.1f} mmol/L\nPente : {pente_str} mmol/L/h ({heure})",
                      priorite="urgent", tags="warning")
                marquer("hyper")

        except KeyboardInterrupt:
            print("\n  Arret du moniteur.")
            notif("Moniteur arrete", "Surveillance glycemie desactivee",
                  priorite="min", tags="stop_sign")
            hc_ping("/fail", f"Arret manuel du moniteur sur {MACHINE}")
            break

        except Exception:
            print("\n  [!] Erreur inattendue :")
            traceback.print_exc()
            hc_ping("/fail", traceback.format_exc()[-2000:])
            try:
                client.login()
                client.fetch_patient()
                print("  [OK] Reconnecte.\n")
            except Exception:
                print("  [X] Reconnexion echouee — nouvel essai au prochain cycle.\n")

        print(f"  (prochaine verification dans {INTERVALLE//60} min)\n")
        time.sleep(INTERVALLE)


if __name__ == "__main__":
    main()