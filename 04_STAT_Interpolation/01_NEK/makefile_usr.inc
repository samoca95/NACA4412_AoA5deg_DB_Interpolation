
INC_SRC := ./inc_src
$(OBJDIR)/misc.o     :./misc.f;                      $(F77) -c $(FL2) -I$(INC_SRC) $< -o $@

