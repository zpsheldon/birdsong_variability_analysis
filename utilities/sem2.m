function y = sem2(x)

y = std(bootstrp(200,@median,x));