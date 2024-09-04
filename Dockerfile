FROM python:3.10 AS base

FROM base AS install
COPY ./src /src
WORKDIR /src
RUN python -m pip install -r requirements.txt

FROM install AS final
COPY ./app /app
WORKDIR /app
ENV DOTNET_SYSTEM_GLOBALIZATION_INVARIANT=1
ENV ASPNETCORE_hostBuilder:reloadConfigOnChange=false
ENV UNITE_COMMAND="python"
ENV UNITE_COMMAND_ARGUMENTS="-u app.py {data}/{proc}"
ENV UNITE_SOURCE_PATH="/src"
ENV UNITE_DATA_PATH="/mnt/data"
ENV UNITE_PROCESS_LIMIT="1"
EXPOSE 80
CMD ["/app/Unite.Commands.Web", "--urls", "http://0.0.0.0:80"]
