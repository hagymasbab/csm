function test_timing()

    c = rand(576);
    v = rand(576,1);    

    a = @() invert_multiply(c,v);
    b = @() invert_multiply2(c,v);
    %b = @() divide(c,vv);
    
    timeit(a)
    timeit(b)
    
end

function invert_multiply(c,v)
    vv = v * v';
    ic = inv(c);
    ic * vv * ic;
end

function invert_multiply2(c,v)
    ic = inv(c);
    m1 = ic * v;
    m1*m1';
end

function divide(c,vv)
    c \ vv / c;
    %eye(size(c,1)) / c;
end