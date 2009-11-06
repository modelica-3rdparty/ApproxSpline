
FC 	= gfortran
FOPTS	= -Wall

LD	= $(FC)

AR	= ar


SRCS	= \
bispev.f  fpader.f   fpcoco.f  fpgrdi.f  fpperi.f  fpsphe.f  polar.f  \
clocur.f  fpadno.f   fpcons.f  fpgrpa.f  fppocu.f  fpsuev.f  profil.f \
cocosp.f  fpadpo.f   fpcosp.f  fpgrre.f  fppogr.f  fpsurf.f  regrid.f \
concon.f  fpback.f   fpcsin.f  fpgrsp.f  fppola.f  fpsysy.f  spalde.f \
concur.f  fpbacp.f   fpcurf.f  fpinst.f  fprank.f  fptrnp.f  spgrid.f \
cualde.f  fpbfout.f  fpcuro.f  fpintb.f  fprati.f  fptrpe.f  sphere.f \
curev.f   fpbisp.f   fpcyt1.f  fpknot.f  fpregr.f  insert.f  splder.f \
curfit.f  fpbspl.f   fpcyt2.f  fpopdi.f  fprota.f  parcur.f  splev.f  \
dblint.f  fpchec.f   fpdeno.f  fpopsp.f  fprppo.f  parder.f  splint.f \
evapol.f  fpched.f   fpdisc.f  fporde.f  fprpsp.f  parsur.f  sproot.f \
          fpchep.f   fpfrno.f  fppara.f  fpseno.f  percur.f  surev.f  \
fourco.f  fpclos.f   fpgivs.f  fppasu.f  fpspgr.f  pogrid.f  surfit.f \

OBJS	=$(SRCS:.f=.o)

all : libdierckx.a

libdierckx.a : $(OBJS)
	ar rvs libdierckx.a $(OBJS)

%.o : %.f
	$(FC) -c $(FOPTS) $<


clean :
	rm -f $(OBJS)

realclean : clean
	rm -f libdierckx.a

