# accrediff/plotting.py
#import cv2 # type: ignore
#from natsort import natsorted # type: ignore
import os
import numpy as np
import matplotlib.pyplot as plt
from .gas_migration import Gas_ModelConfig  # 添加此行
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
def plot_IW_heatmap(
    df_snap_all,
    selected_time,
    Model,
    size_scale   = 2000,
    vmin         = -3,
    vmax         = -0.5,
    xlim         = (0.3, 4),
    ylim         = (0, 0.4),
    sizes_legend = [0.01, 0.1, 1.0],
    figsize      = (16, 10),
    save_path    = None,
):
    """
    绘制指定时刻粒子 ΔIW 热力图。

    Parameters
    ----------
    df_snap_all   : pd.DataFrame  包含 'a', 'e', 'm_e', 'IW' 列
    selected_time : str           时间标签，用于标题显示
    Model         : str           模型名称，用于标题显示
    size_scale    : float         粒子散点大小缩放系数，默认 2000
    vmin / vmax   : float         colorbar 范围（线性刻度），默认 -5 / 1
    xlim / ylim   : tuple         坐标轴范围
    sizes_legend  : list          质量图例节点 (M⊕)
    figsize       : tuple         图像尺寸
    save_path     : str | None    若指定则保存图片，否则直接显示
    """
    # ── 建图 ────────────────────────────────────────────────────
    fig = plt.figure(figsize=figsize)
    gs  = fig.add_gridspec(2, 1, height_ratios=[20, 1], wspace=0.35)
    ax  = fig.add_subplot(gs[0])
    cax = fig.add_subplot(gs[1])

    # ── IW 粒子：线性 Normalize 着色 ────────────────────────────
    scatter = ax.scatter(
        df_snap_all['a'], df_snap_all['e'],
        s          = size_scale * df_snap_all['m_e'],
        c          = df_snap_all['IW'],
        cmap       = 'viridis',
        norm       = plt.Normalize(vmin=vmin, vmax=vmax),
        alpha      = 1.0,
        edgecolors = 'black',
        linewidths = 0.2,
        zorder     = 3,
    )

    # ── Colorbar ─────────────────────────────────────────────────
    cbar = plt.colorbar(scatter, cax=cax, orientation='horizontal', pad=0.0)
    cbar.set_label('bulk oxygen fugacity', fontsize=20, fontweight='bold', labelpad=10)
    cbar.ax.tick_params(labelsize=12)

    # ── 坐标轴格式 ───────────────────────────────────────────────
    ax.set_xlabel('Semi-Major Axis (AU)', fontsize=25, fontweight='bold')
    ax.set_ylabel('Eccentricity',         fontsize=25, fontweight='bold')
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    #ax.grid(True, linestyle='--', alpha=0.3)
    #ax.set_title(f'{Model} | T = {selected_time} | Particles colored by bulk $\\Delta$IW',
    #             fontsize=20, fontweight='bold', pad=20)
    ax.text(
        0.02, 0.98,
        f'T = {selected_time}',
        transform=ax.transAxes,
        fontsize=20, fontweight='bold', va='top', ha='left',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', 
                  alpha=0.7, edgecolor='gray'),
    )
    ax.tick_params(labelsize=13)

    # ── 质量大小图例 ─────────────────────────────────────────────
    size_labels      = [f'{m} M⊕' for m in sizes_legend]
    size_plot_legend = [size_scale * m for m in sizes_legend]
    size_handles = [
        plt.Line2D([0], [0], marker='o', color='w', label=label,
                   markerfacecolor='#555555', markersize=2*np.sqrt(s / np.pi),
                   alpha=0.8, markeredgecolor='black', markeredgewidth=1.2)
        for s, label in zip(size_plot_legend, size_labels)
    ]
    ax.legend(handles=size_handles,
              loc='upper right', fontsize=15, frameon=True,
              title='Particle Mass', title_fontsize=15, edgecolor='black',
              fancybox=True, shadow=True, framealpha=0.95,
              labelspacing=1.4, handleheight=1.4, handletextpad=1.4, borderpad=1.2,
              )

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"✓ Saved: {save_path}")
    else:
        plt.show()
    plt.close(fig)

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