# accrediff/plotting.py
#import cv2 # type: ignore
#from natsort import natsorted # type: ignore
import os
#**************************************************************************************************************************************
def images_to_video(image_folder, output_video, frame_rate):
    try:
        import cv2
        from natsort import natsorted
    except ImportError:
        raise ImportError(
            "OpenCV is required for images_to_video()."
            " Please install it using 'pip install opencv-python'."
        )
    images = [img for img in os.listdir(image_folder) if img.endswith((".png", ".jpg", ".jpeg"))]
    #images.sort()  # Ensure the images are in the correct order
    images = natsorted(images) 

    # Read the first image to get the width and height
    first_image_path = os.path.join(image_folder, images[0])
    first_image = cv2.imread(first_image_path)
    height, width, layers = first_image.shape

    # Define the codec and create VideoWriter object
    fourcc = cv2.VideoWriter_fourcc(*'mp4v')  # Use 'XVID' or 'mp4v' for .mp4
    video = cv2.VideoWriter(output_video, fourcc, frame_rate, (width, height))

    count = 0
    for image in images:
        image_path = os.path.join(image_folder, image)
        img = cv2.imread(image_path)
        if img is not None and img.shape[0] == height and img.shape[1] == width:
            video.write(img)
            count += 1
        else:
            print(f"警告: 读取图片失败或尺寸不一致 {image_path}")

    # Release the video writer
    video.release()
    print(f"视频生成完成，写入帧数: {count}")
#**************************************************************************************************************************************
# 新增：迁移图生成器
def generate_migration_map_v2(
    cfg_template: Gas_ModelConfig,
    a_grid: np.ndarray,
    m_grid: np.ndarray,
    t_code: float = 0.0,
    use_saturation: bool = True
) -> dict:
    """
    生成迁移图：da/dt(a, m_p) 二维数据
    
    Parameters
    ----------
    cfg_template : Gas_ModelConfig
        基础配置模板
    a_grid : np.ndarray
        半长轴网格 [AU]
    m_grid : np.ndarray
        行星质量网格 [M_Earth]
    t_code : float
        代码时间单位
    use_saturation : bool
        是否使用饱和模型
    
    Returns
    -------
    dict with keys:
        - 'a_grid', 'm_grid': 网格坐标
        - 'da_dt': 迁移速率数组 (len(m_grid), len(a_grid))
        - 'color_data': 用于绘图的数据（相对化）
    """
    n_a = len(a_grid)
    n_m = len(m_grid)
    
    da_dt_map = np.zeros((n_m, n_a))
    
    for i, mass in enumerate(m_grid):
        for j, a in enumerate(a_grid):
            # 更新配置
            cfg = Gas_ModelConfig(
                particle_mass_mearth=mass,
                a0=a,
                tau_decay_myr=cfg_template.tau_decay_myr,
                a_alpha1=cfg_template.a_alpha1,
                a_alpha2=cfg_template.a_alpha2,
                Sigma0_cgs=cfg_template.Sigma0_cgs,
            )
            
            # 构建模型并计算 da/dt
            if use_saturation:
                model = Gas_build_model_v2(cfg, use_saturation=True)
            else:
                model = Gas_build_model(cfg)
            
            integrator = model["integrator"]
            da_dt = integrator.da_dt(t_code, a)
            da_dt_map[i, j] = da_dt
    
    # 归一化：用于颜色映射
    # 参考论文：δa/a = da/dt × 1 Myr（归一化到 1 Myr 尺度）
    units = Gas_UnitSystem()
    da_per_myr = da_dt_map * units.code_to_myr(1.0)
    color_data = da_per_myr / np.maximum(np.abs(a_grid[None, :]), 0.01)  # 相对速率
    
    return {
        'a_grid': a_grid,
        'm_grid': m_grid,
        'da_dt': da_dt_map,
        'da_dt_per_myr': da_per_myr,
        'color_data': color_data,
    }


def plot_migration_map_v2(
    migration_map: dict,
    title: str = "Migration Map: da/dt vs (a, m_p)",
    figsize: tuple = (12, 8),
    vmin: float = -3,
    vmax: float = 3,
    cmap: str = "RdBu_r",
    output_path: str = None
) -> None:
    """
    绘制迁移图（类似 Fig. 2）
    
    Parameters
    ----------
    migration_map : dict
        generate_migration_map_v2 的输出
    title : str
        图标题
    figsize : tuple
        图像尺寸
    vmin, vmax : float
        颜色条范围（log 对数，单位：AU/Myr）
    cmap : str
        色彩映射
    output_path : str
        保存路径
    """
    import matplotlib.pyplot as plt
    
    a_grid = migration_map['a_grid']
    m_grid = migration_map['m_grid']
    da_dt_per_myr = migration_map['da_dt_per_myr']
    
    # ---- 对数化处理（为了显示正负速率） ----
    # 方法：sign(da/dt) × log10(|da/dt| + 1e-6)
    da_dt_abs = np.abs(da_dt_per_myr) + 1e-6
    da_dt_signed_log = np.sign(da_dt_per_myr) * np.log10(da_dt_abs)
    
    fig, ax = plt.subplots(figsize=figsize)
    
    # 绘制 2D 数据
    im = ax.pcolormesh(
        a_grid, m_grid, da_dt_signed_log,
        cmap=cmap, vmin=vmin, vmax=vmax,
        shading='auto'
    )
    
    # ---- 添加零迁移线（等高线） ----
    contour = ax.contour(
        a_grid, m_grid, da_dt_per_myr,
        levels=[0], colors='black', linewidths=2
    )
    ax.clabel(contour, inline=True, fontsize=10, fmt='Zero migration')
    
    # ---- 添加快速迁移线（可选） ----
    # 强内向：da/dt = -a（厚实线）
    # 强外向：da/dt = +a（虚线）
    ax.contour(
        a_grid, m_grid, da_dt_per_myr,
        levels=[-a_grid[None, :], a_grid[None, :]],  # 复杂：简化为单一值
        colors='gray', linewidths=1, alpha=0.5
    )
    
    # ---- 轴标签与标题 ----
    ax.set_xlabel('Orbital Radius $a$ (AU)', fontsize=14)
    ax.set_ylabel('Planet Mass (M⊕)', fontsize=14)
    ax.set_title(title, fontsize=15)
    ax.set_yscale('log')
    
    # ---- 颜色条 ----
    cbar = plt.colorbar(im, ax=ax, label=r'$\log_{10}(|da/dt|)$ (AU/Myr)')
    
    plt.tight_layout()
    
    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
    
    plt.show()