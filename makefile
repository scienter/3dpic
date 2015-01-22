EXEC = show
CC = mpicc 
OBJS = main.o parameterSetting.o findparam.o boundary.o loadPlasma3D.o saveFile.o fieldSolve.o loadLaser.o fieldShareY_DSX.o fieldShareZ_DSX.o interpolation.o particlePush.o rearrangeParticles.o particleShareY.o particleShareZ.o removeEdge.o updateCurrent_DSX.o
# loadLaser.o fieldSolve.o fieldShareY_DSX.o interpolation.o particlePush.o updateCurrent_DSX.o rearrangeParticles.o particleShareY.o removeEdge.o movingDomain.o probe.o dumpData.o clean.o filter.o boostShot.o pml.o
INCL = constants.h laser.h mesh.h particle.h plasma.h
LIBS = -lm
$(EXEC):$(OBJS)
	$(CC) $(OBJS) $(LIBS) -o $(EXEC)
$(OBJS):$(INCL)
clean:
	@rm -f *.o *~ $(EXEC)
