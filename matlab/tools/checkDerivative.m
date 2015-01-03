function disc = checkDerivative(formula,deriv,init,builtin)
    
    function [f,g] = fandg(x)
        f = formula(x);
        g = deriv(x);
    end

    if builtin
        
        options = optimoptions('fminunc','GradObj','on','DerivativeCheck','on','FinDiffType','central');
        disc = fminunc(@fandg,init,options);
    
%     expder = deriv(init);
% 
%     function [f,g] = placeholder(x)
%         f = 1;
%         g = zeros(size(expder));
%     end
%     options = optimoptions('fmincon','GradObj','on','DerivativeCheck','on','GradConstr','on');
%     disc = fmincon(@placeholder,init,[],[],[],[],[],[],@fandg,options);
    
    else

        % TODO unsure whether will work in the R^N -> R^M case
      
        initval = formula(init);
        initder = deriv(init);
        outputdim = size(initval(:),1);
        step = 10;
        nPoint = 6;
        disc = [];
        fprintf('Point %d/',nPoint);
        for k=1:nPoint
            printCounter(k);
            step = step*0.1;
            numdiff = zeros([size(init) size(initval)]);
            %size(numdiff)
            for i=1:size(init,1)
                for j=1:size(init,2)
                    act = init;
                    act(i,j) = act(i,j) + step;
                    actval = formula(act);
                    numdiff(i,j,:,:) = (actval-initval)/step;                  
                end                
            end
            summeddiff = numdiff;            
            % summing up differences for output dimensions
            if ndims(summeddiff)>2
                %summeddiff = sum(sum(summeddiff,1),1);            
                summeddiff = squeeze(summeddiff);
            end
%             size(numdiff)
%             size(summeddiff)
%             size(initder)
            disc_mat = (initder-summeddiff).^2;
            discrepancy = sum(disc_mat(:));
            disc=[disc discrepancy];
        end
        fprintf('\n');
        plot(log(disc))
    end
end
            