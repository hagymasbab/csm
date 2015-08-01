function A = gsmLearnA(X,C,x_sigma,z_shape,z_scale,A_init,z_min,z_max,z_res,batch_size,max_step,saveA,randseed)
    setrandseed(randseed);

    z_vals = linspace(z_min,z_max,z_res);
    z_square = z_vals.^2;
    
    z_dens = zeros(batch_size,z_res);
    z_priors = zeros(z_res,1);
    for i = 1:z_res
        z_priors(i,1) = gampdf(z_vals(i),z_shape,z_scale);
    end
    
    iC = inv(C);
    C_u_post = cell(1,z_res);
    mu_u_post = cell(batch_size,z_res);
    
    A = A_init;
    if saveA > 0
        A_iter = cell(1,max_step);
        A_iter{1} = A;
    end
    N_all = size(X,1);
    Du = size(C,1);
    Dx = size(X,2);
    
    skipped_images = 0;
    
    for s = 1:max_step
        %A_prev = A;
        printCounter(s,'maxVal',max_step,'stringVal','EM step');
        
        X_act = X(chooseKfromN(batch_size,N_all),:);
        
        % calculate posterior covariance matrices of U
        for i = 1:z_res
            C_u_post{i} = inv( iC + (z_square(i) / x_sigma^2) * (A'*A) );
        end
        
        % calculate all Z posterior densities and U posterior means for the batch
        for t = 1:batch_size
            xt = X_act(t,:)';
            numerators = zeros(z_res,1);
            for i = 1:z_res
                %like_cov = x_sigma * eye(Dx) + z_square(i) * A*C*A';
                like_cov = x_sigma^2 * eye(Dx) + z_square(i) * A*C*A';
                likelihood = mvnpdf(xt,zeros(size(xt)),like_cov);
                numerators(i,1) = likelihood * z_priors(i,1);
                
%                 % TEST
%                 if isnan(numerators(i,1))
%                     fprintf('e');
%                     numerators(i,1) = 0;
%                 end
                
                mu_u_post{t,i} = (z_vals(i) / x_sigma^2) * C_u_post{i} * A' * xt;
            end
            denominator = sum(numerators);
            
            if denominator == 0              
                z_dens(t,:) = 0;
                skipped_images = skipped_images + 1;
                %fprintf('e');
            else
                for i = 1:z_res
                    z_dens(t,i) = numerators(i,1) / denominator;
                end
            end
        end
        lz = sum(z_dens,1)';                        
        
        % check for NaNs
        
        if any(isnan(z_dens))
            fprintf('z came out nan\n');
        end
        for i=1:z_res
            if any(isnan(C_u_post{i}))
                fprintf('Cu came out nan\n');
            end
            for t=1:batch_size
                if any(isnan(mu_u_post{t,i}))
                    fprintf('mu came out nan\n');
                end
            end
        end
            
        
        % calculate component formulas for A*
        
        psi = zeros(Du);        
        ksi = zeros(Du);
        
        for i = 1:z_res
            psi_part = zeros(Du);
            ksi_part = zeros(Du);
            for t = 1:batch_size
                xt = X_act(t,:)';
                mu_act = mu_u_post{t,i};
                p_mu = z_dens(t,i) * mu_act;
                psi_part = psi_part + xt * p_mu';
                ksi_part = ksi_part + p_mu * mu_act';
            end
            psi = psi + z_vals(i) * psi_part;
            ksi = ksi + z_square(i) * (C_u_post{i} * lz(i) + ksi_part);
        end        
        
        A = psi / ksi;
        
        if rem(s,saveA)==0
            A_iter{s+1} = A;
            save('gsm_A_iter.mat','A_iter','C','x_sigma','z_shape','z_scale','-v7.3');
        end
    end
    fprintf('Skipped %d images because posterior of z came out flat for them',skipped_images);
end

