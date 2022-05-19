import open3d as o3d

pcd = o3d.io.read_point_cloud("zivid_bot_view_1.ply")
o3d.visualization.draw_geometries([pcd])
