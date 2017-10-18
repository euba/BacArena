library(BacArena)

set.seed(5000)
data(Ec_core)

bac = Bac(model=Ec_core)
arena = Arena(n=20, m=2)
arena <- addOrg(arena, bac, amount=1)
arena <- addDefaultMed(arena, bac)

sim <- simEnv(arena, time = 2)
sim@exchanges
