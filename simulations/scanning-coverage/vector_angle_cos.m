function cos_theta = vector_angle_cos(u, v)
cos_theta = dot(u, v) / sqrt(sum(u.^2) * sum(v.^2));
end