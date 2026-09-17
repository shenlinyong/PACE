FROM python:3.12-slim

ARG PACE_EXTRAS=io,ml
ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

WORKDIR /opt/pace
COPY . .
RUN case ",${PACE_EXTRAS}," in \
        *,sequence,*) python -m pip install --no-cache-dir 'torch>=2.6,<3' --index-url https://download.pytorch.org/whl/cpu ;; \
    esac \
    && python -m pip install --no-cache-dir ".[${PACE_EXTRAS}]" \
    && useradd --system --uid 10001 --create-home pace \
    && mkdir -p /work \
    && chown pace:pace /work

USER pace
WORKDIR /work
ENTRYPOINT ["PACE"]
CMD ["--help"]
