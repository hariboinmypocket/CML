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



# sys/stat.h: "file is dataless object". A file iCloud has evicted keeps every
# byte of its metadata, so st_size, mtime and os.path.exists all still answer
# correctly and open() succeeds -- only the first read blocks, waiting on a
# download. With the disk full that download cannot land, so the read never
# returns and never raises.
SF_DATALESS = 0x40000000


class EvictedFileError(OSError):
    """An input file's contents live in iCloud, not on this disk."""


def _file_flags(path: Path) -> int | None:
    """st_flags for `path`, or None where the platform has no such field."""
    try:
        return os.stat(path).st_flags
    except AttributeError:  # pragma: no cover - st_flags is BSD/macOS only
        return None


def require_local(path: str | Path) -> Path:
    """Refuse to read a file whose contents iCloud has evicted.

    Reading one does not fail, it hangs forever, which is the least visible
    failure this project has hit. It cost a session to diagnose: the mirror's
    flags read "compressed,dataless" with zero allocated blocks while `ls`
    still reported all 109,830,144 bytes, so a stalled read looked like failing
    hardware rather than a file that simply was not here. `du` reporting the
    file as 0B was the real signature and got dismissed as a glitch.

    Checked by flag, not by st_blocks: APFS compression also reports zero
    blocks, and most of this repo's source is compressed and perfectly local.
    """
    target = Path(path)
    flags = _file_flags(target)
    if flags is not None and flags & SF_DATALESS:
        raise EvictedFileError(
            f"{target} is in iCloud, not on this disk (dataless), so reading "
            f"it would hang rather than fail. Recover it with:\n"
            f"    brctl download '{target}'\n"
            f"Turn off iCloud Drive's \"Optimize Mac Storage\" to stop macOS "
            f"evicting this project's inputs when the disk fills."
        )
    return target
