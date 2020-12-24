function e = wientropy(s)

e = mean(log(max(s,.0001))) - log(mean(max(s,.0001)));