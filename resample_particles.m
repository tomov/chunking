function new_particles = resample_particles(particles, w)

    w=w/sum(w);
    new_particles = particles;
    for j=1:length(particles)
        k=find(mnrnd(1,w));
        new_particles(j)=particles(k);
    end
