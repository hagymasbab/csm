function test_timing()

    c = rand(576);
    vv = c * c';

    a = @() invert_multiply(c,vv);
    b = @() divide(c,vv);
    
    timeit(a)
    timeit(b)
    
end

function invert_multiply(c,vv)
    ic = inv(c);
    ic * vv * ic;
end

function divide(c,vv)
    c \ vv / c;
    %eye(size(c,1)) / c;
end