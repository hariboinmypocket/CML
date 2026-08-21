from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path
from typing import Iterator

import mysql.connector
from mysql.connector import MySQLConnection

from .config import Settings


def connect(settings: Settings, include_database: bool = True) -> MySQLConnection:
    options = {
        "host": settings.host,
        "port": settings.port,
        "user": settings.user,
        "password": settings.password,
        "autocommit": False,
        "charset": "utf8mb4",
        "use_unicode": True,
    }
    if include_database:
        options["database"] = settings.database
    return mysql.connector.connect(**options)


def initialize_database(settings: Settings, schema_path: str | Path) -> None:
    with connect(settings, include_database=False) as connection:
        cursor = connection.cursor()
        safe_name = settings.database.replace("`", "``")
        cursor.execute(
            f"CREATE DATABASE IF NOT EXISTS `{safe_name}` "
            "CHARACTER SET utf8mb4 COLLATE utf8mb4_0900_ai_ci"
        )
        cursor.execute(f"USE `{safe_name}`")
        sql = Path(schema_path).read_text(encoding="utf-8")
        cursor.execute(sql, map_results=True)
        while cursor.nextset():
            pass
        connection.commit()
        cursor.close()


@contextmanager
def transaction(settings: Settings) -> Iterator[MySQLConnection]:
    connection = connect(settings)
    try:
        yield connection
        connection.commit()
    except Exception:
        connection.rollback()
        raise
    finally:
        connection.close()
