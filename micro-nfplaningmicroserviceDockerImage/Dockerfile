FROM matlabruntime/r2023a/release/update6/308000000000000000

COPY ./applicationFilesForMATLABCompiler /usr/bin/mlrtapp
RUN chmod -R a+rX /usr/bin/mlrtapp/*

RUN useradd -ms /bin/bash appuser
USER appuser

ENTRYPOINT ["/opt/matlabruntime/R2023a/bin/glnxa64/muserve", "-a", "/usr/bin/mlrtapp/NFPSarchive.ctf"]
