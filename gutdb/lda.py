"""A fixed, bounded representation of LDA scores for association features.

Softsign x/(1+abs(x)) is monotonic, symmetric, and independent of other rows.
It preserves comparison-order sign; it is not a probability, a fold change,
or a correction for between-study differences. Rounding can create ties.
"""
from __future__ import annotations

from decimal import Decimal, ROUND_HALF_UP


COLUMN = "lda_score_normalized"
EXPRESSION = """CASE WHEN effect_type = 'LDA' THEN
    ROUND(effect_size / (1.0 + ABS(effect_size)), 3)
    ELSE NULL END"""
SUMMARY_SELECT = """SELECT t.id AS taxon_id, t.scientific_name,
    d.id AS disease_id, d.name AS disease_name, d.mesh_id,
    a.direction, a.effect_type, a.effect_size, s.project_accession,
    a.lda_score_normalized
FROM taxon_disease_associations a
JOIN taxa t ON t.id = a.taxon_id
JOIN diseases d ON d.id = a.disease_id
JOIN phenotype_comparisons c ON c.id = a.comparison_id
JOIN studies s ON s.id = c.study_id"""


def normalize_lda_score(value: float | None) -> Decimal | None:
    """Reference calculation, with symmetric rounding to three decimals."""
    if value is None:
        return None
    x = Decimal(str(value))
    if not x.is_finite():
        raise ValueError("LDA score must be finite")
    return (x / (1 + abs(x))).quantize(Decimal("0.001"), rounding=ROUND_HALF_UP)


def install_normalized_lda(connection, *, sqlite: bool = False) -> bool:
    """Add a generated column and expose it in the summary view.

    Existing effect sizes are never updated. Generated columns automatically
    follow later inserts/updates without changing any ingestion code. MySQL
    DDL commits implicitly; SQLite callers control their transaction.
    Returns whether a column was added. Safe to repeat with this definition.
    """
    cursor = connection.cursor()
    try:
        if sqlite:
            cursor.execute("PRAGMA table_xinfo(taxon_disease_associations)")
            exists = any(row[1] == COLUMN for row in cursor.fetchall())
            number_type = "REAL"
        else:
            cursor.execute("""SELECT COLUMN_NAME, EXTRA
                FROM information_schema.COLUMNS
                WHERE TABLE_SCHEMA = DATABASE()
                  AND TABLE_NAME = 'taxon_disease_associations'
                  AND COLUMN_NAME = %s""", (COLUMN,))
            row = cursor.fetchone()
            exists = row is not None
            if exists and "GENERATED" not in row[1]:
                raise ValueError(f"Existing {COLUMN} is not generated; inspect before migration")
            number_type = "DECIMAL(4,3)"
        if not exists:
            cursor.execute(f"""ALTER TABLE taxon_disease_associations
                ADD COLUMN {COLUMN} {number_type}
                GENERATED ALWAYS AS ({EXPRESSION}) VIRTUAL""")
        if sqlite:
            cursor.execute("DROP VIEW IF EXISTS v_taxon_disease_summary")
            cursor.execute("CREATE VIEW v_taxon_disease_summary AS " + SUMMARY_SELECT)
        else:
            cursor.execute("CREATE OR REPLACE VIEW v_taxon_disease_summary AS " + SUMMARY_SELECT)
        return not exists
    finally:
        cursor.close()
