import os
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo, ZoneInfoNotFoundError


def _zoneinfo_from_name(name):
    if not name:
        return None

    zone_name = name.strip()
    if zone_name.startswith(":"):
        zone_name = zone_name[1:]
    if not zone_name or zone_name.startswith("/"):
        return None

    try:
        return ZoneInfo(zone_name)
    except ZoneInfoNotFoundError:
        return None


def _zoneinfo_from_localtime_symlink():
    localtime = Path("/etc/localtime")
    try:
        target = localtime.resolve()
    except OSError:
        return None

    parts = target.parts
    if "zoneinfo" not in parts:
        return None

    zoneinfo_index = parts.index("zoneinfo")
    zone_name = "/".join(parts[zoneinfo_index + 1 :])
    if zone_name.startswith(("posix/", "right/")):
        zone_name = "/".join(zone_name.split("/")[1:])
    return _zoneinfo_from_name(zone_name)


def get_local_timezone():
    for env_var in ("AIDAMRI_TIMEZONE", "TZ"):
        timezone = _zoneinfo_from_name(os.environ.get(env_var))
        if timezone is not None:
            return timezone

    for timezone_file in (Path("/etc/timezone"), Path("/var/db/zoneinfo")):
        try:
            timezone = _zoneinfo_from_name(timezone_file.read_text().strip())
        except OSError:
            timezone = None
        if timezone is not None:
            return timezone

    timezone = _zoneinfo_from_localtime_symlink()
    if timezone is not None:
        return timezone

    return datetime.now().astimezone().tzinfo

