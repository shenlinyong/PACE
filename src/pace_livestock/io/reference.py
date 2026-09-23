"""Reference FASTA access for experimental annotation and CpG validation."""

from pathlib import Path

from ..errors import PaceError


class Reference:
    """Read regular uncompressed FASTA by byte offsets without loading the genome."""

    def __init__(self, path):
        self.path = Path(path)
        self.index = {}
        self.handle = open(path, "rb")
        chrom, offset, length, bases, width, short = None, 0, 0, 0, 0, False
        try:
            while True:
                position = self.handle.tell()
                line = self.handle.readline()
                if not line or line.startswith(b">"):
                    if chrom is not None:
                        if not length:
                            raise PaceError(f"Empty FASTA sequence: {chrom}")
                        self.index[chrom] = (offset, length, bases, width)
                    if not line:
                        break
                    chrom = line[1:].split()[0].decode()
                    if chrom in self.index:
                        raise PaceError(f"Duplicate FASTA chromosome: {chrom}")
                    offset, length, bases, width, short = self.handle.tell(), 0, 0, 0, False
                else:
                    sequence = line.rstrip(b"\r\n")
                    if chrom is None or not sequence or set(sequence.upper()) - set(b"ACGTN"):
                        raise PaceError("Reference must be uncompressed FASTA containing A/C/G/T/N")
                    if not bases:
                        bases, width, offset = len(sequence), len(line), position
                    elif (
                        short
                        or len(sequence) > bases
                        or (len(sequence) == bases and len(line) != width and line.endswith(b"\n"))
                    ):
                        raise PaceError("FASTA lines must have regular width except the final line")
                    short = len(sequence) < bases or not line.endswith(b"\n")
                    length += len(sequence)
            if not self.index:
                raise PaceError("Empty FASTA")
        except Exception:
            self.handle.close()
            raise

    @property
    def sizes(self):
        return {k: v[1] for k, v in self.index.items()}

    def fetch(self, chrom, start, end):
        if chrom not in self.index or start < 0 or end < start or end > self.index[chrom][1]:
            raise PaceError(f"FASTA interval outside reference: {chrom}:{start}-{end}")
        offset, _, bases, width = self.index[chrom]
        result = []
        while start < end:
            n = min(end - start, bases - start % bases)
            self.handle.seek(offset + (start // bases) * width + start % bases)
            result.append(self.handle.read(n).decode("ascii").upper())
            start += n
        return "".join(result)

    def close(self):
        self.handle.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
