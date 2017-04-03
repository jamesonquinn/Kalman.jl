function ppath(p)
  if LOAD_PATH[end] != p
    push!(LOAD_PATH,p)
  else
    println(p," already in LOAD_PATH.")

  end
end

ppath("/Users/chema/mydev/Kalman.jl/src")
ppath("/Users/chema/mydev/")
ppath("/Users/chema/mydev/Gadfly.jl/src")

macro load(pkg)
  quote
    if isdefined(Symbol($(string(pkg))))
      #try
      #  reload(ucfirst(string(pkg)))
      #catch
        reload($(string(pkg)))
      #end
    else
      import $pkg
    end
  end
end

@load Gadfly



@load bkf

x0 = bkf.State([1.,1],[1. -.99; -.99 1])
a = [1. .5; -.5 1]
a = a / det(a)
g = [0.05 0; 0 .01]
q = eye(2)
f = bkf.LinearModel(a,g,q)
h = eye(2)
g = [0.15 0; 0 .01]
z = bkf.LinearObservationModel(h)
kf0 = bkf.BasicKalmanFilter(x0,f,z)

pf = bkf.toParticleSet(kf0,1000)
pf1 = bkf.toParticleSet(kf0,1)

dt = .1
t = collect(0:dt:1)

truth = Vector{bkf.ParticleSet}(0)#length(t))
observations = Vector{bkf.Observation}(0)#length(t))
kfs = Vector{bkf.BasicKalmanFilter}(0)#length(t))
pfs = Vector{bkf.ParticleStep}(0)#length(t))
m = [cos(tt) for tt in t]

kf = kf0
push!(kfs, kf)
ps = bkf.ParticleStep(pf)
push!(pfs, ps)
pf

push!(truth, pf1)
push!(observations, bkf.Observation(pf1,1))









kf1 = bkf.predictupdate(kf,observations[1])
bkf.predictupdate!(kf,observations[1])
#@test kf1.x == kf.x

for i in 2:length(t)-1
    push!(truth, bkf.ap(truth[i-1]))
    push!(observations, bkf.Observation(truth[i],1))
    kf2 = bkf.predict(kf)
    kf = bkf.update(kf2,observations[i])
    push!(kfs, kf)
    ps = bkf.ParticleStep(ps,observations[i])
    push!(pfs, ps)
end



using DataFrames

nshow = 30
df1 = DataFrame(id=[x for x in pfs[end].r.s], time=length(t),
          x=vec(pfs[end].p.particles[1,:]), y=vec(pfs[end].p.particles[2,:]))
#using Gadfly
df2 = DataFrame(id=[x for x in pfs[end].r.s], time=length(t)- 1,
          x=vec(pfs[end-1].p.particles[1,pfs[end].r.s]),
          y=vec(pfs[end-1].p.particles[2,pfs[end].r.s]))

df = vcat(df1[1:nshow,:],df2[1:nshow,:])
dfh = hcat(df1,df2)

pp=Gadfly.plot(df,x="x",y="y",group="id",Gadfly.Geom.line,
    Gadfly.Guide.annotation(
        Gadfly.compose(Gadfly.context(),
        Gadfly.ellipse(-.25,-0.5,.05,.10), Gadfly.fill(nothing),
        Gadfly.stroke("red"))) )
inch = Gadfly.inch
pp
show(pp)
Gadfly.draw(Gadfly.SVG("/Users/chema/mydev/myplot2.svg", 4inch, 3inch), pp)

pp2=Gadfly.plot(dfh,x="x",y="x_1",Gadfly.Geom.point)
Gadfly.draw(Gadfly.SVG("/Users/chema/mydev/myplot3.svg", 4inch, 3inch), pp2)

3









#room for atom error display









@load bkf
ds = [10,25,50,100]
ln = size(ns,1)
reps = 5
finkelmean = zeros(ln,reps)
frankenmean = zeros(ln,reps)
partmean = zeros(ln,reps)
idealmean = zeros(ln,reps)
for width in 1:ln

    d = ds[width]
    bleed = .25
    jitter = .1
    jitterbleed = .1 # ends up being like twice this, because hits on left and right, blech.
    temper = .9
    basenoise = 1
    lownoise = .04
    lowgap = 3
    noisebleed = .1
    mciter = 6

    
    x0 = bkf.State(zeros(d),eye(d))

    a = SymTridiagonal(ones(d),bleed * ones(d-1))
    a = a / det(a)^(1/d)

    b = SymTridiagonal(ones(d),-jitterbleed * ones(d-1))
    b[1,1] = b[d,d] = 1 - jitterbleed^2 / (1-jitterbleed^2)
    g = inv(b)
    g = g * (temper / det(g))^(1/d)

    q = eye(d)

    f = bkf.LinearModel(a,g,q)


    h = eye(d) * basenoise
    for i in 1:lowgap:d
        h[i,i] = lownoise
    end

    b = SymTridiagonal(ones(d),-noisebleed * ones(d-1))
    b[1,1] = b[d,d] = 1 - noisebleed^2 / (1-noisebleed^2)
    r = inv(b)
    r = r / det(r)^(1/d)

    z = bkf.LinearObservationModel(h)
    kf0 = bkf.BasicKalmanFilter(x0,f,z)


    nfp,npf,nfapf = (5, 25, 25)
    #nfp,npf,nfapf = (31, 1000, 200)
    nfp,npf,nfapf = (100, 10000, 2000)

    for r in 1:reps


        fpf = bkf.toParticleSet(kf0,nfp)
        pf = bkf.toParticleSet(kf0,npf)
        fapf = bkf.toFrankenSet(kf0,n,5)
        fp = bkf.FinkelToe(fpf)

        pf1 = bkf.toParticleSet(kf0,1)

        T = 3
        t = collect(0:T)

        truth = Vector{bkf.ParticleSet}(0)#length(t))
        observations = Vector{bkf.Observation}(0)#length(t))
        kfs = Vector{bkf.BasicKalmanFilter}(0)#length(t))
        pfs = Vector{bkf.ParticleStep}(0)#length(t))
        fps = Vector{bkf.AbstractFinkel}(0)#length(t))
        faps = Vector{bkf.FrankenStep}(0)

        kf = kf0
        push!(kfs, kf)
        ps = bkf.ParticleStep(pf)
        push!(pfs, ps)
        push!(fps, fp)
        fap = bkf.FrankenStep(fapf)
        push!(faps, fap)

        push!(truth, pf1)
        push!(observations, bkf.Observation(pf1,1))



        for i in 2:length(t)-1
            push!(truth, bkf.ap(truth[i-1]))
            push!(observations, bkf.Observation(truth[i],1))
            kf2 = bkf.predict(kf)
            kf = bkf.update(kf2,observations[i])
            push!(kfs, kf)
            ps = bkf.ParticleStep(ps,observations[i])
            push!(pfs, ps)
            fp = bkf.FinkelParticles(fp, observations[i], 4)
            push!(fps, fp)
            fap = bkf.FrankenStep(fap, observations[i])
            push!(faps, fap)
        end

        dif = kfs[end].x.p - kfs[end].x.p'
        mean(dif)
        mean(kfs[end].x.p)
        finalDist = bkf.toDistribution(kfs[end])
        finkelmean[width,r] = mean(log(bkf.pdf(finalDist,fps[end].tip.particles)))/d
        partmean[width,r] = mean(log(bkf.pdf(finalDist,pfs[end].p.particles)))/d
        frankenmean[width,r] = mean(log(bkf.pdf(finalDist,faps[end].p.particles)))/d
        idealmean[width,r] = mean(log(bkf.pdf(finalDist,rand(finalDist,50))))/d
    end
end

[idealmean,finkelmean,frankenmean,partmean]










#....
