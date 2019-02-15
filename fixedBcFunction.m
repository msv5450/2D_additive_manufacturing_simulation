function f = fixedBcFunction(x,y,value)

nn = length(x);
f = ones(nn,1) * value;