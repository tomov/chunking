function results = forward(T, num_particles, init_fn, choice_fn, update_fn)

    choices = [];
    for i=1:num_particles
        particles(i) = init_fn();
        [w(i) c(i)] = choice_fn(1, particles(i));
    end
    liks(1) = mean(w);
    choices = [choices; mean(c)]; % TOTE this only works for Bernoulli choices!

    particles = resample_particles(particles, w);

    for t=2:T
        for i=1:num_particles
            particles(i) = update_fn(t-1, particles(i));
            [w(i) c(i)] = choice_fn(t, particles(i));
        end
        liks(n) = mean(w);
        choices = [choices; mean(c)]; % TOTE this only works for Bernoulli choices!
    
        particles = resample_particles(particles, w);
    end
    
    for i=1:num_particles
        particles(i) = update_fn(T, particles(i)); % last update of posterior
    end

    results.liks = liks;

    %save forward.mat
