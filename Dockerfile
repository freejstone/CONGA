FROM python:3.9
LABEL maintainer="Michael Riffle <mriffle@uw.edu>"
ENV APPDIR=/app

WORKDIR  $APPDIR
ADD entrypoint.sh /usr/local/bin/entrypoint.sh
ADD . .

RUN chmod 755 /usr/local/bin/entrypoint.sh && \
    pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir .

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD []
