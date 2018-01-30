AR = ar -ru
RM = rm -f

CFLAGS = -O3
CPPFLAGS = -DHAVE_UNISTD_H=1 -DHAVE_STDARG_H=1 -DHAVE_HIDDEN=1 -DHAVE_MEMCPY=1 -DHAVE_MEMMOVE=1 -DNDEBUG
INC = -I. -I"dierckx" -I"modelica"

TARGETDIR = linux64

OBJS = \
	dierckx/bispev.o \
	dierckx/clocur.o \
	dierckx/cocosp.o \
	dierckx/concon.o \
	dierckx/concur.o \
	dierckx/cualde.o \
	dierckx/curev.o \
	dierckx/curfit.o \
	dierckx/dblint.o \
	dierckx/evapol.o \
	dierckx/fourco.o \
	dierckx/fpader.o \
	dierckx/fpadno.o \
	dierckx/fpadpo.o \
	dierckx/fpback.o \
	dierckx/fpbacp.o \
	dierckx/fpbfout.o \
	dierckx/fpbisp.o \
	dierckx/fpbspl.o \
	dierckx/fpchec.o \
	dierckx/fpched.o \
	dierckx/fpchep.o \
	dierckx/fpclos.o \
	dierckx/fpcoco.o \
	dierckx/fpcons.o \
	dierckx/fpcosp.o \
	dierckx/fpcsin.o \
	dierckx/fpcurf.o \
	dierckx/fpcyt1.o \
	dierckx/fpcyt2.o \
	dierckx/fpdeno.o \
	dierckx/fpdisc.o \
	dierckx/fpfrno.o \
	dierckx/fpgivs.o \
	dierckx/fpgrdi.o \
	dierckx/fpgrpa.o \
	dierckx/fpgrre.o \
	dierckx/fpgrsp.o \
	dierckx/fpinst.o \
	dierckx/fpintb.o \
	dierckx/fpknot.o \
	dierckx/fpopdi.o \
	dierckx/fpopsp.o \
	dierckx/fporde.o \
	dierckx/fppara.o \
	dierckx/fppasu.o \
	dierckx/fpperi.o \
	dierckx/fppocu.o \
	dierckx/fppogr.o \
	dierckx/fppola.o \
	dierckx/fprank.o \
	dierckx/fprati.o \
	dierckx/fpregr.o \
	dierckx/fprota.o \
	dierckx/fprppo.o \
	dierckx/fprpsp.o \
	dierckx/fpseno.o \
	dierckx/fpspgr.o \
	dierckx/fpsphe.o \
	dierckx/fpsuev.o \
	dierckx/fpsurf.o \
	dierckx/fpsysy.o \
	dierckx/fptrnp.o \
	dierckx/fptrpe.o \
	dierckx/insert.o \
	dierckx/parcur.o \
	dierckx/parder.o \
	dierckx/parsur.o \
	dierckx/percur.o \
	dierckx/pogrid.o \
	dierckx/polar.o \
	dierckx/profil.o \
	dierckx/regrid.o \
	dierckx/spalde.o \
	dierckx/spgrid.o \
	dierckx/sphere.o \
	dierckx/splder.o \
	dierckx/splev.o \
	dierckx/splint.o \
	dierckx/surev.o \
	dierckx/surfit.o \
	mdc_dierckx.o

ALL_OBJS = $(OBJS)

all: clean libs

libs: libmdc_dierckx.a
	cp $^ ../Library/$(TARGETDIR)

libmdc_dierckx.a: $(OBJS)
	$(AR) $@ $(OBJS)

%.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) $(INC) -c -o $@ $<

clean:
	$(RM) $(ALL_OBJS)
	$(RM) *.a
	$(RM) ../Library/$(TARGETDIR)/*.a
	$(RM) ../Library/$(TARGETDIR)/$(TARGETDIR).tar.xz