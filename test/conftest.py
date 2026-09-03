"""
Make mrd2recon importable for the tests, with or without the mrd package installed.

mrd2recon imports mrd at module scope, but the arithmetic these tests exercise never
touches it: only one mrd attribute is evaluated at import time on python 3.14 (emit's
`array_type=mrd.ArrayType.USER_MAP` default argument), plus the annotation names on 3.12,
where annotations are still evaluated eagerly. So a small stand-in is enough to reach every
numeric function, and the suite runs on numpy and scipy alone.

The real package is preferred whenever it is installed, so this never masks a genuine
change in the mrd schema. To run against it, use an interpreter that has it, e.g.
    ~/.local/share/mamba/envs/mrd/bin/python -m unittest discover -s test -t .

Imported for its side effects by the test modules; also picked up automatically as a
conftest if pytest is ever added.
"""

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


class _Enum:
    """Stands in for an mrd enum: any attribute is a distinct, printable sentinel."""

    def __init__(self, name):
        self._name = name
        self._members = {}

    def __getattr__(self, item):
        if item.startswith("_"):
            raise AttributeError(item)
        return self._members.setdefault(item, f"{self._name}.{item}")


class _Record:
    """Stands in for an mrd record type: keyword arguments become attributes."""

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


def _install_mrd_stub():
    import types

    mrd = types.ModuleType("mrd")
    for name in ("ArrayType", "ArrayDimension", "ArrayImageType", "AcquisitionFlags"):
        setattr(mrd, name, _Enum(name))

    class _StreamItem:
        NdArrayDouble = staticmethod(lambda a: ("NdArrayDouble", a))
        NdArrayComplexDouble = staticmethod(lambda a: ("NdArrayComplexDouble", a))
        Acquisition = staticmethod(lambda a: ("Acquisition", a))

    mrd.StreamItem = _StreamItem

    class _ArrayMetaValue:
        String = staticmethod(lambda v: ("s", v))
        Int64 = staticmethod(lambda v: ("i", v))
        Float64 = staticmethod(lambda v: ("f", v))

    mrd.ArrayMetaValue = _ArrayMetaValue
    mrd.ArrayMeta = dict

    for name in ("NdArray", "NdArrayHeader", "Header", "Acquisition",
                 "UserParametersType", "UserParameterDoubleType", "UserParameterLongType",
                 "BinaryMrdReader", "BinaryMrdWriter"):
        setattr(mrd, name, type(name, (_Record,), {}))

    sys.modules["mrd"] = mrd
    return mrd


try:  # pragma: no cover - depends on the environment, both branches are exercised in CI-less use
    import mrd  # noqa: F401
    USING_REAL_MRD = True
except ImportError:
    _install_mrd_stub()
    USING_REAL_MRD = False
