%This function uses the implemented python function "face_of_span" by the
%liberary cdd inside matlab
function U = face_of_span(V)
    pyrun("import numpy as np");
    V_np = pyrun("VGI_np = np.array(vec)", "VGI_np", vec = V);
    pyrun("from span_to_face import face_of_span");
    U = double(pyrun("U = face_of_span(V)", "U", V = V_np));
end
