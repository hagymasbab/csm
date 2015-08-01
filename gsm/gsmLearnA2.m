function A = gsmLearnA2(X,C,x_sigma,z_shape,z_scale,A_init,z_res,batch_size,max_step,saveA,randseed)
    setrandseed(randseed);    
    N_all = size(X,1);
    Du = size(C,1);
    Dx = size(X,2);
    if size(A_init,1) ~= Dx
        error('Data is of dimension %d and filterbank output dimension is %d',Dx,size(A_init,1));
    end
    if size(A_init,2) ~= Du
        error('Latent covariance dimension %d and filterbank intput dimension is %d',Du,size(A_init,2));
    end
    
    iC = inv(C);        
    noiseCov = eye(Dx) * x_sigma^2;
    x0 = zeros(Dx,1);
    
    A = A_init;
    if saveA > 0
        A_iter = {};
        A_iter{1} = A;
    end    
            
    x_indices = [];
    
    for s = 1:max_step        
        printCounter(s,'maxVal',max_step,'stringVal','EM step');
        
        act_x_indices = chooseKfromN(batch_size,N_all);
        x_indices = [x_indices; act_x_indices];
        save('xind.mat','x_indices');
        
        X_act = X(act_x_indices,:);        
        sATA = (A' * A) / x_sigma^2;
        ACAT = A * C * A';     
        
        % csalas
        ACAT = nearestSPD(ACAT);
        
        psi = zeros(Du);
        ksi = zeros(Dx,Du);
        for t = 1: batch_size
            %printCounter(t,'maxVal',batch_size,'stringVal','Observation');
            xt = X_act(t,:)';
            
            % get a MAP estimate of z
            [z_map,z_min,z_max] = get_pz_x_max(xt, noiseCov, ACAT, x0, z_shape, z_scale);
            z_vals = linspace(z_min,z_max,z_res);
            z_square = z_vals.^2;
            
            % calculate posterior density of z
            z_post_numerators = zeros(z_res,1);
            for i = 1:z_res
                like_cov = noiseCov + z_square(i) * ACAT;
                likelihood = mvnpdf(xt,zeros(size(xt)),like_cov);
                prior = gampdf(z_vals(i),z_shape,z_scale);
                z_post_numerators(i,1) = likelihood * prior;
            end
            z_post_denominator = sum(z_post_numerators);
            
            psi_part = zeros(Du,1);            
            
            for i = 1:z_res
                C_u_post = inv( iC + z_square(i) * sATA );
                mu_u_post = (z_vals(i) / x_sigma^2) * C_u_post * A' * xt;
                z_post_dens = z_post_numerators(i) / z_post_denominator;
                
                psi_part = psi_part + z_post_dens * z_vals(i) * mu_u_post;
                ksi = ksi + z_square(i) * z_post_dens * ( C_u_post + mu_u_post * mu_u_post' );
            end
            psi = psi + xt * psi_part';
        end
        
        A = psi / ksi;
        
        if rem(s,saveA)==0
            A_iter{end+1} = A;
            save('gsm_A_iter.mat','A_iter','C','x_sigma','z_shape','z_scale','-v7.3');
        end
        
    end
    
end

