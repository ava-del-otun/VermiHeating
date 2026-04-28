from __future__ import annotations

import argparse
import datetime as dt
import json
import statistics
import urllib.parse
import urllib.request
from pathlib import Path


DEFAULT_MONTHLY_MEAN_C = [-1.34, 0.08, 3.89, 8.91, 14.78, 19.74, 23.32, 21.98, 18.16, 12.44, 5.54, 0.53]
DEFAULT_EVENT_SOURCE_TEMPLATE = (
    "https://www.ncei.noaa.gov/access/services/data/v1"
    "?dataset=daily-summaries"
    "&stations={station}"
    "&startDate={start_date}"
    "&endDate={end_date}"
    "&dataTypes=TMAX,TMIN,TAVG"
    "&units=standard"
    "&format=json"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Derive monthly heat-wave spell anomaly bounds from NOAA daily-summaries data. "
            "Spell anomalies are computed at spell level (one anomaly per spell), grouped by spell start month."
        )
    )
    parser.add_argument("--station", default="USW00014740", help="NOAA station id (default: USW00014740)")
    parser.add_argument("--start-date", default="2016-01-01", help="Start date (YYYY-MM-DD)")
    parser.add_argument("--end-date", default="2025-12-31", help="End date (YYYY-MM-DD)")
    parser.add_argument(
        "--event-source-url",
        default="",
        help="Optional explicit NOAA daily-summaries URL. If omitted, one is built from station/start/end.",
    )
    parser.add_argument(
        "--input-json",
        default="",
        help="Optional local NOAA daily-summaries JSON file to use instead of HTTP fetch.",
    )
    parser.add_argument(
        "--save-source-json",
        default="",
        help="Optional path to save the NOAA daily-summaries rows used for derivation.",
    )
    parser.add_argument(
        "--output-json",
        default="",
        help="Optional path to write derived monthly spell-anomaly statistics as JSON.",
    )
    parser.add_argument(
        "--monthly-mean-c",
        default=",".join(str(v) for v in DEFAULT_MONTHLY_MEAN_C),
        help="Comma-separated 12-value monthly baseline mean-temperature array in C.",
    )
    parser.add_argument("--min-spell-length-days", type=int, default=3, help="Minimum consecutive hot days per spell.")
    parser.add_argument("--tmax-threshold-f", type=float, default=90.0, help="Heat-wave threshold on daily TMAX (F).")
    return parser.parse_args()


def parse_float(value) -> float | None:
    if value is None:
        return None
    if isinstance(value, (float, int)):
        return float(value)
    text = str(value).strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def fahrenheit_to_celsius(value_f: float) -> float:
    return (float(value_f) - 32.0) * 5.0 / 9.0


def parse_monthly_mean_array(text: str) -> list[float]:
    parts = [segment.strip() for segment in str(text).split(",")]
    values = [float(part) for part in parts if part]
    if len(values) != 12:
        raise ValueError(f"--monthly-mean-c must contain 12 numeric values; received {len(values)} values.")
    return values


def build_event_source_url(station: str, start_date: str, end_date: str) -> str:
    return DEFAULT_EVENT_SOURCE_TEMPLATE.format(
        station=urllib.parse.quote_plus(str(station).strip()),
        start_date=urllib.parse.quote_plus(str(start_date).strip()),
        end_date=urllib.parse.quote_plus(str(end_date).strip()),
    )


def load_source_rows_from_json(path: Path) -> list[dict]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if isinstance(payload, list):
        return [dict(row) for row in payload if isinstance(row, dict)]
    if isinstance(payload, dict):
        for key in ("rows", "data", "results", "items"):
            value = payload.get(key)
            if isinstance(value, list):
                return [dict(row) for row in value if isinstance(row, dict)]
    raise ValueError(f"Unsupported NOAA source JSON shape in {path}")


def fetch_source_rows_from_url(url: str) -> list[dict]:
    with urllib.request.urlopen(url, timeout=120) as response:
        payload = json.loads(response.read().decode("utf-8"))
    if isinstance(payload, list):
        return [dict(row) for row in payload if isinstance(row, dict)]
    if isinstance(payload, dict):
        for key in ("rows", "data", "results", "items"):
            value = payload.get(key)
            if isinstance(value, list):
                return [dict(row) for row in value if isinstance(row, dict)]
    raise ValueError("Unsupported NOAA response payload shape.")


def parse_daily_records(rows: list[dict]) -> list[dict]:
    records: list[dict] = []
    for row in rows:
        date_raw = row.get("DATE", row.get("date", ""))
        date_text = str(date_raw).strip()
        if not date_text:
            continue
        try:
            day = dt.date.fromisoformat(date_text[:10])
        except ValueError:
            continue
        tmax_f = parse_float(row.get("TMAX", row.get("tmax")))
        tmin_f = parse_float(row.get("TMIN", row.get("tmin")))
        tavg_f = parse_float(row.get("TAVG", row.get("tavg")))
        if tavg_f is None and tmax_f is not None and tmin_f is not None:
            tavg_f = 0.5 * (tmax_f + tmin_f)
        records.append(
            {
                "date": day,
                "tmax_f": tmax_f,
                "tmin_f": tmin_f,
                "tavg_f": tavg_f,
            }
        )
    records.sort(key=lambda item: item["date"])
    return records


def detect_heatwave_spells(records: list[dict], min_spell_length_days: int, tmax_threshold_f: float) -> list[list[dict]]:
    spells: list[list[dict]] = []
    run: list[dict] = []
    minimum_days = max(int(min_spell_length_days), 1)
    threshold_f = float(tmax_threshold_f)

    for row in records:
        tmax_f = row.get("tmax_f")
        is_hot_day = tmax_f is not None and float(tmax_f) >= threshold_f
        if is_hot_day:
            if run and (row["date"] - run[-1]["date"]).days != 1:
                if len(run) >= minimum_days:
                    spells.append(run)
                run = []
            run.append(row)
        else:
            if len(run) >= minimum_days:
                spells.append(run)
            run = []
    if len(run) >= minimum_days:
        spells.append(run)
    return spells


def round_monthly(values: list[float], digits: int = 2) -> list[float]:
    return [float(round(float(v), digits)) for v in values]


def derive_monthly_spell_anomaly_bounds(
    records: list[dict],
    monthly_mean_c: list[float],
    min_spell_length_days: int,
    tmax_threshold_f: float,
) -> dict:
    spells = detect_heatwave_spells(records, min_spell_length_days=min_spell_length_days, tmax_threshold_f=tmax_threshold_f)
    spell_anomalies_by_start_month: list[list[float]] = [[] for _ in range(12)]

    for spell in spells:
        if not spell:
            continue
        start_month_idx = int(spell[0]["date"].month) - 1
        day_anomalies: list[float] = []
        for day in spell:
            tavg_f = day.get("tavg_f")
            if tavg_f is None:
                continue
            day_month_idx = int(day["date"].month) - 1
            monthly_ref_c = float(monthly_mean_c[day_month_idx])
            anomaly_c = max(fahrenheit_to_celsius(float(tavg_f)) - monthly_ref_c, 0.0)
            day_anomalies.append(float(anomaly_c))
        if day_anomalies:
            spell_anomaly_c = float(sum(day_anomalies) / len(day_anomalies))
            spell_anomalies_by_start_month[start_month_idx].append(spell_anomaly_c)

    mean_vals: list[float] = []
    min_vals: list[float] = []
    max_vals: list[float] = []
    std_vals: list[float] = []
    spell_counts: list[int] = []

    for month_vals in spell_anomalies_by_start_month:
        spell_counts.append(int(len(month_vals)))
        if month_vals:
            mean_vals.append(float(sum(month_vals) / len(month_vals)))
            min_vals.append(float(min(month_vals)))
            max_vals.append(float(max(month_vals)))
            std_vals.append(float(statistics.pstdev(month_vals)) if len(month_vals) > 1 else 0.0)
        else:
            mean_vals.append(0.0)
            min_vals.append(0.0)
            max_vals.append(0.0)
            std_vals.append(0.0)

    return {
        "spellCountTotal": int(sum(spell_counts)),
        "spellCountByStartMonth": spell_counts,
        "meanSpellAnomalyAboveMonthlyMean_C": round_monthly(mean_vals),
        "minSpellAnomalyAboveMonthlyMean_C": round_monthly(min_vals),
        "maxSpellAnomalyAboveMonthlyMean_C": round_monthly(max_vals),
        "spellAnomalyStd_C": round_monthly(std_vals),
    }


def main() -> None:
    args = parse_args()
    monthly_mean_c = parse_monthly_mean_array(args.monthly_mean_c)
    source_rows: list[dict]
    resolved_source = ""
    source_mode = ""

    if str(args.input_json).strip():
        source_path = Path(str(args.input_json).strip()).expanduser().resolve()
        source_rows = load_source_rows_from_json(source_path)
        resolved_source = str(source_path)
        source_mode = "local_json"
    else:
        source_url = str(args.event_source_url).strip() or build_event_source_url(
            station=args.station,
            start_date=args.start_date,
            end_date=args.end_date,
        )
        source_rows = fetch_source_rows_from_url(source_url)
        resolved_source = source_url
        source_mode = "http_json"

    records = parse_daily_records(source_rows)
    derived = derive_monthly_spell_anomaly_bounds(
        records,
        monthly_mean_c=monthly_mean_c,
        min_spell_length_days=int(args.min_spell_length_days),
        tmax_threshold_f=float(args.tmax_threshold_f),
    )
    derived.update(
        {
            "station": str(args.station),
            "startDate": str(args.start_date),
            "endDate": str(args.end_date),
            "tmaxThreshold_F": float(args.tmax_threshold_f),
            "minSpellLength_days": int(args.min_spell_length_days),
            "monthlyMean_C": [float(v) for v in monthly_mean_c],
            "dailyMeanProxy": "TAVG if present else 0.5*(TMAX+TMIN) from NOAA daily-summaries",
            "spellGroupingBasis": "group by spell start month; per-spell anomaly is mean of daily anomalies across spell days",
            "sourceMode": source_mode,
            "sourceResolved": resolved_source,
            "recordCount": int(len(records)),
            "generatedUtc": dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z"),
        }
    )

    if str(args.save_source_json).strip():
        save_path = Path(str(args.save_source_json).strip()).expanduser().resolve()
        save_path.parent.mkdir(parents=True, exist_ok=True)
        save_path.write_text(json.dumps(source_rows, indent=2), encoding="utf-8")
        print(f"Saved NOAA rows: {save_path}")

    if str(args.output_json).strip():
        output_path = Path(str(args.output_json).strip()).expanduser().resolve()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(json.dumps(derived, indent=2), encoding="utf-8")
        print(f"Saved derived bounds: {output_path}")

    print(f"Source mode: {source_mode}")
    print(f"Source: {resolved_source}")
    print(f"Daily records parsed: {len(records)}")
    print(f"Detected heat-wave spells: {derived['spellCountTotal']}")
    print("Derived monthly spell anomaly arrays (C):")
    print(json.dumps(derived, indent=2))


if __name__ == "__main__":
    main()
