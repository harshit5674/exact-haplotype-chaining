"""Memory-light independent graph-major verifier (gap_validation_light.cpp; V*H distance words instead of V*V).

Same interface and checks as gap_validation.py; not part of the solver or its timing."""
from array import array
import ctypes as ct
import hashlib
from pathlib import Path
import subprocess
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
_LIB=None

def buf(a):
    t={'i':ct.c_int,'q':ct.c_int64}[a.typecode]
    return (t*len(a)).from_buffer(a) if a else None

def build():
    global _LIB
    if _LIB is not None:return _LIB
    src=HERE/'gap_validation_light.cpp';flags=['-std=c++20','-O2','-shared','-fPIC']
    target=ROOT/'results/native_build'/('gap-validation-light-'+hashlib.sha256(src.read_bytes()+repr(flags).encode()).hexdigest()[:16]+('.dylib' if sys.platform=='darwin' else '.so'))
    if not target.exists():subprocess.run(['clang++',*flags,str(src),'-o',str(target)],check=True)
    lib=ct.CDLL(str(target));ip=ct.POINTER(ct.c_int);qp=ct.POINTER(ct.c_int64);vp=ct.c_void_p
    lib.validation_index.argtypes=[ct.c_int,ct.c_int,ct.c_int,qp,ip,ip,ip];lib.validation_index.restype=vp
    lib.validation_batch.argtypes=[vp,ct.c_int,qp,ct.c_int64];lib.validation_batch.restype=vp
    lib.validation_check.argtypes=[vp,qp,qp,qp,qp];lib.validation_check.restype=ct.c_int
    lib.validation_index_free.argtypes=[vp];lib.validation_batch_free.argtypes=[vp]
    lib.validation_error.restype=ct.c_char_p
    _LIB=lib;return lib

class Validator:
    def __init__(self,graph):
        self.lib=build();paths=graph['paths'];flat=array('i');offset=array('i',[0])
        for p in paths:flat.extend(p);offset.append(len(flat))
        lengths=array('q',graph['lengths']);edges=array('i',[x for e in graph['edges'] for x in e])
        self.handle=self.lib.validation_index(len(lengths),len(paths),len(edges)//2,buf(lengths),buf(edges),buf(offset),buf(flat))
        if not self.handle:raise AssertionError(self.lib.validation_error().decode())
    def batch(self,anchors,gamma):return ValidationBatch(self,anchors,gamma)
    def close(self):
        if self.handle:self.lib.validation_index_free(self.handle);self.handle=None
    def __del__(self):self.close()

class ValidationBatch:
    def __init__(self,index,anchors,gamma):
        self.index=index;flat=array('q',[x for a in anchors for x in (a.v,a.x,a.y,a.c,a.d,a.weight)])
        self.handle=index.lib.validation_batch(index.handle,len(anchors),buf(flat),gamma)
        if not self.handle:raise AssertionError(index.lib.validation_error().decode())
    def check(self,result):
        if self.index.lib.validation_check(self.handle,buf(result.C),buf(result.B),buf(result.R),buf(result.back)):
            raise AssertionError(self.index.lib.validation_error().decode())
    def close(self):
        if self.handle:self.index.lib.validation_batch_free(self.handle);self.handle=None
    def __del__(self):self.close()
