function y = Clip(X, lowerValue, upperValue)
    y = min(max(X, lowerValue), upperValue);
end