%%
function X_double = datetime2double(X)
    X_double = fix(X.Year*1e8 + X.Month*1e6 + X.Day*1e4 + X.Hour*1e2  + X.Minute);
end