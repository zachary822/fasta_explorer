from pydantic import BaseModel, AfterValidator
from typing import Annotated
import re

DNA_REGEX = re.compile(r"^[=ACMGRSVTWYHKDBN]*$")


def validate_dna(code: str):
    if not DNA_REGEX.match(code):
        raise ValueError("invalid DNA sequence")

    return code


class Sequence(BaseModel):
    sequence: Annotated[str, AfterValidator(validate_dna)]
    reverse_complement: str
    gc_fraction: float


class FileResult(BaseModel):
    sequences: list[Sequence]
