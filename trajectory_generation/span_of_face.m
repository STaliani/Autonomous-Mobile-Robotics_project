%This function uses the implemented python function "span_of_face" by the
%liberary cdd inside matlab
function V = span_of_face(U)
    pyrun("import numpy as np");
    U_np = pyrun("U_np = np.array(vec)", "U_np", vec = U);
    pyrun("from face_to_span import span_of_face");
    V = double(pyrun("V = span_of_face(U)", "V", U = U_np));
end