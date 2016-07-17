# Top level makefile
install:    build distclean

build:
	@echo "Building software dependence..."
	@	if test -f makedeps.sh; \
		then \
			sh makedeps.sh; \
			else echo "Fail to build dependence!"; \
			exit; \
		fi
	$(MAKE) -C src build

clean:
	@echo "Cleaning software..."
	$(MAKE) -C src clean

distclean:
	@echo "Removing installation files..."
	$(MAKE) -C src distclean
	-/bin/rm -rf config.* configure makedeps.sh Makefile autom4te.cache
