function disc = checkDerivative(formula,deriv,init)
    initval = formula(init);
    initder = deriv(init);
    step = 0.1;
    disc = [];
    for k=1:6
        step = step*0.1;
        discrepancy = 0;
        for i=1:size(init,1)
            for j=1:size(init,2)
                act = init;
                act(i,j) = act(i,j) + step;
                actval = formula(act);
                numdiff = (actval-initval)/step;
                discrepancy = discrepancy + initder(i,j)-numdiff;
            end
        end
        disc=[disc discrepancy];
    end
    plot(disc)
end
            