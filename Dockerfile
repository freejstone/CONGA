FROM python:3.9
LABEL maintainer="Michael Riffle <mriffle@uw.edu>"

ADD entrypoint.sh /usr/local/bin/entrypoint.sh
ADD CONGA.py /usr/local/bin/CONGA.py
ADD peptides.py /usr/local/bin/peptides.py
ADD requirements.txt /usr/local/bin/requirements.txt

RUN chmod 755 /usr/local/bin/entrypoint.sh

RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r /usr/local/bin/requirements.txt

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD []
