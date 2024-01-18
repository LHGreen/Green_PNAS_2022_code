function Rsq = Rsquare(y, yCalc)

Rsq = 1-sum((y-yCalc).^2)/sum((y-mean(y)).^2);

end