from pydantic import BaseModel


class Sequence(BaseModel):
    sequence: str
    reverse_complement: str
    gc_fraction: float


class FileResult(BaseModel):
    sequences: list[Sequence]
