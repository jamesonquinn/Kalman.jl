
truth[end].particles[10:15,1] - [mean(rsamps[i,:]) for i = 10:15]
truth[end].particles[10:15,1] - [mean(pfs[end].p.particles[i,:]) for i = 10:15]
truth[end].particles[10:15,1] - [mean(fps[end].tip.particles[i,:]) for i = 10:15]
truth[end].particles[10:15,1] - [mean(faps[end].p.particles[i,:]) for i = 10:15]
ruth[2].filter.z.h[10:15,10:15]
cov(faps[2].p.particles[10:15,:],2)
cov(rsamps[10:15,:],2)






gf = Gadfly
gf.plot(gf.layer(x=rsamps[1,1:100],y=rsamps[2,1:100],
                gf.Geom.point,gf.Theme(default_color=gf.color("red"))),
       gf.layer(x=fra1.base[1,1:100],y=fra1.base[2,1:100],
                gf.Geom.point,gf.Theme(default_color=gf.color("blue"))),
       gf.layer(x=fra1.tip.particles[1,1:100],y=fra1.tip.particles[2,1:100],
                gf.Geom.point,gf.Theme(default_color=gf.color("green"))),
       gf.layer(x=fra1k.tip.particles[1,1:100],y=fra100.tip.particles[2,1:100],
                gf.Geom.point,gf.Theme(default_color=gf.color("purple")))
    )

mean(fra1.base,2)


mean(fra1.tip.particles,2)


mean(fra100.tip.particles,2)


mean(fra1k.tip.particles,2)


mean(rsamps,2)




cov(fra1.tip.particles,2)


cov(fra1k.tip.particles,2)


cov(rsamps,2)



fra1k.base
