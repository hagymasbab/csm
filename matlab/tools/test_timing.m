function test_timing()

    c = rand(576);
    c = c*c';
    v = rand(576,1);    

    a = @() invert_multiply(c,v);
    b = @() invert_multiply2(c,v);
    %b = @() divide(c,vv);
    
    r1 = invert_multiply(c,v);
    r2 = invert_multiply2(c,v);
    isequal(r1,r2)
    viewImage(r1-r2,'useMax',false);
    
    timeit(a)
    timeit(b)
    
end

function ret = invert_multiply(c,v)
    vv = v * v';
    ic = inv(c);
    ret = ic * vv * ic;
end

function ret = invert_multiply2(c,v)
    ic = inv(c);
    m1 = ic * v;
    ret = m1*m1';
end

function divide(c,vv)
    c \ vv / c;
    %eye(size(c,1)) / c;
end