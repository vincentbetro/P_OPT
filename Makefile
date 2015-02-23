
SRCDIRS = SMOOTH MAIN

all:
	@for dir in $(SRCDIRS); do \
	  cd $$dir ; \
	  pwd ; \
	  $(MAKE) ; \
	  cd .. ; \
	done

clean:
	@for dir in $(SRCDIRS); do \
	  cd $$dir ; \
	  pwd ; \
	  $(MAKE) clean ; \
	  cd .. ; \
	done
	/bin/rm -f P_OPT.*

cleanest:
	@for dir in $(SRCDIRS); do \
	  cd $$dir ; \
	  pwd ; \
	  $(MAKE) clean ; \
	  cd .. ; \
	done
	/bin/rm -f P_OPT.*
	/bin/rm -f P_OPT_out.*
	/bin/rm -f *.jou
	/bin/rm -f *.uns
	/bin/rm -f *.params
