#!/bin/sh

rm -f config.cache acconfig.h

if [ ! -d config ]; then
		echo "mkdir config"
		mkdir config
fi
autoreconf --force --install -I config && exit 0

exit 1
