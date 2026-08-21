from __future__ import annotations

from dataclasses import dataclass
import os
from pathlib import Path


def load_env(path: str | Path = ".env") -> None:
    """Load simple KEY=VALUE settings without requiring python-dotenv."""
    env_path = Path(path)
    if not env_path.exists():
        return
    for raw_line in env_path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        key = key.strip()
        value = value.strip().strip('"').strip("'")
        os.environ.setdefault(key, value)


@dataclass(frozen=True)
class Settings:
    host: str
    port: int
    user: str
    password: str
    database: str
    batch_size: int
    gmrepo_base_url: str

    @classmethod
    def from_env(cls, env_file: str | Path = ".env") -> "Settings":
        load_env(env_file)
        return cls(
            host=os.getenv("GUTDB_HOST", "127.0.0.1"),
            port=int(os.getenv("GUTDB_PORT", "3306")),
            user=os.getenv("GUTDB_USER", "root"),
            password=os.getenv("GUTDB_PASSWORD", ""),
            database=os.getenv("GUTDB_NAME", "gut_microbiome"),
            batch_size=int(os.getenv("GUTDB_BATCH_SIZE", "500")),
            gmrepo_base_url=os.getenv(
                "GMREPO_BASE_URL", "https://gmrepo.humangut.info"
            ).rstrip("/"),
        )

