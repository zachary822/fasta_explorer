import io
from pathlib import Path

import bam_reader
import bam_reader.cutils
import bam_reader.utils
import Bio.SeqIO
from fastapi import FastAPI, HTTPException, UploadFile
from models import Sequence, FileResult

app = FastAPI()


@app.post("/upload")
async def upload(file: UploadFile) -> FileResult:
    match Path(file.filename).suffix:
        case ".fasta":
            sequences = []

            for seq in Bio.SeqIO.parse(io.TextIOWrapper(file.file), "fasta"):
                sequences.append(
                    Sequence(
                        sequence=str(seq.seq),
                        reverse_complement=str(seq.seq.reverse_complement()),
                        gc_fraction=bam_reader.cutils.gc_fraction(bytes(seq)),
                    )
                )

            return FileResult(sequences=sequences)
        case ".bam":
            sequences = []

            buff = io.BytesIO(bam_reader.decompress_bzgf(file.file))

            for seq in bam_reader.extract_sequence(buff):
                sequences.append(
                    Sequence(
                        sequence=seq.decode("ascii"),
                        reverse_complement=bam_reader.utils.reverse_complement(
                            seq
                        ).decode("ascii"),
                        gc_fraction=bam_reader.cutils.gc_fraction(seq),
                    )
                )
            return FileResult(sequences=sequences)

    raise HTTPException(status_code=400, detail="Invalid file format")
