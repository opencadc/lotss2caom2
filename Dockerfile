FROM python:3.10-slim

RUN apt-get update --no-install-recommends && \
    apt-get install -y \
        build-essential \
        curl \
        less \
        libaio-dev \
        unzip

WORKDIR /usr/src/app

RUN mkdir -p /usr/local/src
# ADD https://lta.lofar.eu/software/instantclient-basic-linux.x64-12.2.0.1.0.zip /usr/local/src
ADD https://lta.lofar.eu/software/instantclient-basic-linux-11.2.0.3.0.zip /usr/local/src
RUN cd /usr/local/src && \
    unzip instantclient-basic-linux-11.2.0.3.0.zip

ADD https://lta.lofar.eu/software/instantclient-sdk-linux.x64-11.2.0.3.0.zip /usr/local/src
RUN cd /usr/local/src && \
    unzip instantclient-sdk-linux.x64-11.2.0.3.0.zip

ADD https://lta.lofar.eu/software/lofar_lta-2.8.0.tar.gz /usr/local/src
RUN cd /usr/local/src && \
    tar xvf lofar_lta-2.8.0.tar.gz 

RUN ln -s /usr/local/src/instantclient_11_2/libclntsh.so.11.1 /usr/local/src/instantclient_11_2/libclntsho.so
ENV ORACLE_HOME /usr/local/src/instantclient_11_2
ENV LD_LIBRARY_PATH /usr/local/lib/instantclient_11_2:/usr/local/src/instantclient_11_2
RUN echo "/usr/local/src/instantclient_11_2" > /etc/ld.so.conf.d/oracle.conf && ldconfig
RUN mkdir /usr/local/src/instantclient_11_2/log

RUN cd /usr/local/src/lofar_lta-2.8.0 && \
    python setup.py install_oracle && \
    python setup.py install

ARG FINDER_PATH=/usr/local/lib/python${OPENCADC_PYTHON_VERSION}/site-packages/footprintfinder.py
ADD http://www.eso.org/~fstoehr/footprintfinder.py ${FINDER_PATH}
RUN chmod 644 ${FINDER_PATH}

ARG OPENCADC_BRANCH=main
ARG OPENCADC_REPO=opencadc

RUN git clone https://github.com/${OPENCADC_REPO}/caom2tools.git && \
    cd caom2tools && \
    git checkout ${OPENCADC_BRANCH} && \
    pip install ./caom2utils && \
    cd ..

RUN pip install git+https://github.com/${OPENCADC_REPO}/caom2pipe@${OPENCADC_BRANCH}#egg=caom2pipe

RUN pip install git+https://github.com/${OPENCADC_REPO}/lotss2caom2@${OPENCADC_BRANCH}#egg=lotss2caom2

COPY ./docker-entrypoint.sh /usr/local/bin

RUN useradd --create-home --shell /bin/bash cadcops
RUN chown -R cadcops:cadcops /usr/src/app
USER cadcops

ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]
